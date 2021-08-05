##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################

"""
===========================
pipeline_QAPA
===========================



:Author: Jack Riley
:Release: $Id$
:Date: |today|
:Tags: Python



Overview
========

This pipeline functions to quantify levels of alternative polyadenylation (APA)
using the QAPA package released by Ha, Blencowe and Morris (2018). The output
from QAPA consists of a database containing individual sample APA quantifications.
Custom python and R scripts are used to interogate this database and subsequently
allow for multiple comparisons to be made based on the design.tsv (see below). 


Usage
=====

To execute this pipeline in an interactive ShARC session or from a personal
device (not recommended) use the following command:

"python <path_to_pipeline_folder>/pipeline_QAPA.py make full -v5"

To execute this pipeline to the ShARC cluster, use the following:

"submit_pipeline <path_to_pipeline_folder>/pipeline_QAPA.py make full -v5"


Configuration
=============

The final stages of pipeline_QAPA require a design.tsv in order to associate file names
to experimental conditions, and subsequently group and compare the outputs of QAPA
between conditions.


Input files
===========

Raw RNA sequencing reads are given as input in .fastq format. Both single-end and paired-end
reads are accepted by the pipeline. Reference transcriptome and design.tsv are also required.


Requirements
============

Packages required:
    - cgat/cgatcore
    - samtools 
    - stringtie
    - portcullis
    - salmon >v1.0 Selective Alignment update
    - qapa =v1.3.1 (or newer if still compatible) --> https://github.com/morrislab/qapa

These packages are installed in the conda environment "qapa-env".
R packages for final analysis/reports are installed in "qapa-env-R". 


Pipeline output
===============

Raw QAPA database will be stored in the ./QAPA/output/raw/pau_results.txt. Grouped and processed
results will be output as an R-markdown document in ./QAPA/output/final_report.Rmd with supporting
datasets saved in this subdirectory in .csv or .tsv format. 


Code
====

"""

###################
##### imports #####
###################

from ruffus import *
from ruffus.combinatorics import product
import sys
import os
import shutil
from gffutils import DataIterator as DataIterator
import sqlite3
import subprocess
import glob
import csv
from cgatcore import experiment as E
import cgat.Sra as Sra
from cgatcore import pipeline as P
import cgatpipelines.tasks.rnaseq as RnaSeq
import tempfile

############################
#####       main       #####
############################

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    'genesets',
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

PARAMS["project_src"]=os.path.dirname(__file__)

RnaSeq.PARAMS = PARAMS


########################
##### QAPA scripts #####
########################
    
@follows(mkdir("QAPA", "QAPA/prereqs"))
@originate("QAPA/prereqs/ensembl_identifiers.txt")
def downloadEnsemblMetadata(outfile):
    '''Download Ensembl gene metadata table via MySQL'''

    ensembl_mart_assembly_version = PARAMS["qapa_ensembl_mart_assembly_version"]
    ensembl_mart_species = PARAMS["qapa_ensembl_mart_species"]

    statement = ''' mysql --user=anonymous 
                        --host=martdb.ensembl.org 
                        --port=5316
                        -A %(ensembl_mart_assembly_version)s
                        -e "select stable_id_1023 as 'Gene stable ID',
                        stable_id_1066 as 'Transcript stable ID',
                        biotype_1020 as 'Gene type',
                        biotype_1064 as 'Transcript type',
                        display_label_1074 as 'Gene name'
                        from %(ensembl_mart_species)s"
                        > %(outfile)s;
                    '''

    P.run(statement, job_condaenv="qapa-env")


@follows(mkdir("QAPA", "QAPA/prereqs"))
@originate("QAPA/prereqs/gencode.basic.txt")
def downloadGencodeAnnotation(outfile):
    '''Download GENCODE gene prediction annotation table via MySQL'''

    gencode_gpat = PARAMS["qapa_gencode_pred_table"]
    gencode_assembly = PARAMS["qapa_gencode_assembly"]

    statement = ''' mysql --user=genome 
                        --host=genome-mysql.cse.ucsc.edu 
                        -A
                        -e "select * from %(gencode_gpat)s" 
                        %(gencode_assembly)s
                        > %(outfile)s;
                    '''

    P.run(statement, job_condaenv="qapa-env")


@follows(mkdir("QAPA", "QAPA/prereqs"))
@originate("QAPA/prereqs/pas_db.bed")
def downloadPASdb(outfile):
    '''Download PolyA site database'''

    pas_db_href = PARAMS["qapa_pas_db_href"]
    pas_db_name = PARAMS["qapa_pas_db_name"]
    unzipped_pas_db_name = pas_db_name[:-3]

    statement = ''' wget %(pas_db_href)s;
                    gzip -d %(pas_db_name)s;
                    mv %(unzipped_pas_db_name)s %(outfile)s;
                    '''

    P.run(statement, job_condaenv="qapa-env")


@follows(mkdir("QAPA", "QAPA/prereqs"))
@originate("QAPA/prereqs/gencode.polyA_sites.bed")
def downloadGencodePolyA(outfile):
    '''Download GENCODE polyA sites track'''

    gencode_polyA = PARAMS["qapa_gencode_polyA_track"]
    gencode_assembly = PARAMS["qapa_gencode_assembly"]

    statement = ''' mysql --user=genome 
                        --host=genome-mysql.cse.ucsc.edu 
                        -A
                        -e "select chrom, txStart, txEnd,
                        name2, score, strand from
                        %(gencode_polyA)s where name2 = 'polyA_site'"
                        -N %(gencode_assembly)s
                        > %(outfile)s;
                    '''

    P.run(statement, job_condaenv="qapa-env")


@follows(downloadEnsemblMetadata, downloadGencodeAnnotation, downloadPASdb, downloadGencodePolyA)
def downloadQAPAprereqs():
    pass


@follows(downloadQAPAprereqs)
@merge([downloadEnsemblMetadata,
        downloadGencodePolyA,
        downloadPASdb,
        downloadGencodeAnnotation],
        "QAPA/utr_lib.bed")
def build3UTRlib(infiles, outfile):
    '''Build 3'UTR library from prerequisites'''

    ensembl_metadata, gencode_polyA, PASdb, gencode_annotation = infiles

    job_threads=2
    job_memory="32G"

    statement = ''' qapa build --db %(ensembl_metadata)s
                        -g %(gencode_polyA)s
                        -p %(PASdb)s
                        %(gencode_annotation)s >
                        %(outfile)s
                    '''

    P.run(statement, job_condaenv="qapa-env", job_threads=job_threads, job_memory=job_memory)

@follows(build3UTRlib)
@transform(build3UTRlib, filter=suffix("lib.bed"), output="seqs.fa")
def extract3UTRseq(infile, outfile):
    '''Extract 3'UTR sequences from reference genome'''
    
    job_threads=4
    job_memory="32G"

    ref_genome_dir = PARAMS["genome_dir"]
    ref_genome = PARAMS["genome"]

    statement = ''' qapa fasta -f %(ref_genome_dir)s/%(ref_genome)s
                        %(infile)s
                        %(outfile)s
                    '''
    
    P.run(statement, job_condaenv="qapa-env", job_memory=job_memory,)

@follows(extract3UTRseq)
@transform(extract3UTRseq, filter=suffix("_seqs.fa"), output=".salmon.index")
def makeSalmonIndex(infile, outfile):
    '''Create salmon index against different 3'UTR isoforms'''

    job_threads = 2
    job_memory = "64G"

    statement = ''' salmon index -t %(infile)s
                        -i %(outfile)s
                    '''

    P.run(statement, job_condaenv="qapa-env", job_memory=job_memory, job_threads=job_threads)

##########################
##### quantification #####
##########################

@follows(mkdir("QAPA/quantification"), mkdir("sorted_bams"),
         makeSalmonIndex)
@transform("input_assemble.dir/*.bam", regex("(.+)/(.+).bam"), output=r"QAPA/quantification/\2/quant.sf")
def quantifyWithSalmon(infile, outfile):
    '''Quantify existing samples against genesets'''
    job_threads=2
    job_memory="32G"

    basefile = os.path.basename(infile)
    sample_name = basefile.split(os.extsep, 1)
    salmonIndex = "QAPA/utr.salmon.index"

    sorted_bam="sorted_bams/" + sample_name[0] + "_sorted.bam"
    fastq1 = P.snip(outfile, "/quant.sf")+".1.fastq"
    fastq2 = P.snip(outfile, "/quant.sf")+".2.fastq"
    fastq0 = P.snip(outfile, "/quant.sf")+".0.fastq"
    outfile = P.snip(outfile, "/quant.sf")

    statement = ''' samtools sort -n %(infile)s -o %(sorted_bam)s;
                    samtools fastq
                        -1 %(fastq1)s
                        -2 %(fastq2)s
                        -0 %(fastq0)s -s /dev/null -n -F 0x900
                        %(sorted_bam)s;
                    paired_end_reads=$(samtools view -c -f 1 %(sorted_bam)s);
                    if [ $paired_end_reads = 0 ]; then
                        salmon quant -i %(salmonIndex)s
                            --libType IU
                            -r %(fastq0)s
                            -o %(outfile)s;
                    else
                        salmon quant -i %(salmonIndex)s
                            --libType IU
                            -1 %(fastq1)s
                            -2 %(fastq2)s
                            -o %(outfile)s;
                    fi; 
                    rm %(fastq1)s; rm %(fastq2)s; rm %(fastq0)s; rm %(sorted_bam)s 
                    '''
    P.run(statement, job_condaenv="qapa-env", job_memory=job_memory, job_threads=job_threads)

@follows(quantifyWithSalmon, mkdir("QAPA/outputs"))
def quant3UTRusage():
    '''Combines quantifications and computes relative proportions for each'''

    statement = "qapa quant --db QAPA/prereqs/ensembl_identifiers.txt QAPA/quantification/*/quant.sf > QAPA/outputs/pau_results.txt"
    os.system(statement)

#################################################
# Comparing QAPA outputs based on design matrix #
#################################################

@follows(quant3UTRusage)
@transform("design.tsv", suffix("design.tsv"), "QAPA/outputs/comparisons_to_make.txt") 
def analyseDesignMatrix(infile, outfile):
    infile = open(infile)
    design = csv.reader(infile, delimiter="\t")
    outfile = open(outfile, "w")
    current_file = __file__ 
    pipeline_path = os.path.abspath(current_file)
    pipeline_directory = os.path.dirname(pipeline_path)
    script_path = "pipeline_QAPA/"
    Rmd_path = os.path.join(pipeline_directory, script_path)
    to_cluster=False


    for row in design:
        variables = len(row) - 1 
        if(variables == 1):
            folder_name = str(row[0] + "_vs_" + str(row[1]))
            info = "~var1"
            outfile.write(folder_name + "\n")
            statement = (""" mkdir QAPA/outputs/""" + folder_name + """ && cp """ + Rmd_path + 
                         """compare_QAPA_outputs.Rmd QAPA/outputs/""" + folder_name + """/compare_""" + 
                         folder_name + """.Rmd""")
            os.system(statement)

        else:
            raise valueError("design.tsv is not configured correctly, see pipeline_QAPA/example_design.tsv for setup")

    infile.close
    outfile.close()

@follows(analyseDesignMatrix)
@transform("QAPA/outputs/*/compare_*.Rmd", suffix(".Rmd"), ".html")
def run_compareQAPA_script(infile, outfile):
    '''Compares outputs of QAPA based on design.tsv input'''
    job_threads = 2
    job_memory = "16G"

    statement = "Rscript -e 'rmarkdown::render(\"%(infile)s\")'"
    P.run(statement, job_condaenv="R-rstudio", job_memory=job_memory)


###################
##### utility #####
###################

@follows(downloadQAPAprereqs, build3UTRlib, extract3UTRseq, makeSalmonIndex, 
        quantifyWithSalmon, quant3UTRusage, analyseDesignMatrix,
        run_compareQAPA_script)
def full():
    pass

##################
###### misc ######
##################

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))