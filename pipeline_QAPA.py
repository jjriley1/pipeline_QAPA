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

##########################################################
#####       main - substantial amount of this        ##### 
#####  code comes from pipeline_utrons/pipeline_DTU  #####
##########################################################

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

###################
##### utility #####
###################

def connect():
    '''utility function to connect to database.
    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.
    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

STRINGTIE_QUANT_FILES=["i_data.ctab", "e_data.ctab", "t_data.ctab",
                       "i2t.ctab", "e2t.ctab"]

######################################
##### assemble novel transcripts #####
######################################

@follows(mkdir("assembled_transcripts.dir"), mkdir("portcullis"))
@transform(["input_assemble.dir/*.bam",
            "input_assemble.dir/*.remote"],
           formatter(),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"])),
           "assembled_transcripts.dir/{basename[0]}.gtf.gz")
def assembleWithStringTie(infiles, outfile):

    infile, reference = infiles
    basefile = os.path.basename(infile)
    job_threads = PARAMS["stringtie_threads"]
    job_memory = PARAMS["stringtie_memory"]
    tmpfile = P.get_temp_filename()
    if os.path.exists(tmpfile): 
    	os.unlink(tmpfile)
    
    statement =  '''portcullis full 
			            -t 1
                        -o portcullis/%(basefile)s/
                        -r %(portcullis_bedref)s
                        -b
                        %(portcullis_fastaref)s
                        %(infile)s &&
                    mv portcullis/%(basefile)s/portcullis.filtered.bam %(tmpfile)s &&
                    rm -r portcullis/%(basefile)s/ &&
                    stringtie %(tmpfile)s
                        -p %(stringtie_threads)s
                        -G <(zcat %(reference)s)
                        %(stringtie_options)s
                        2> %(outfile)s.log
                    | gzip > %(outfile)s &&
                    rm %(tmpfile)s'''

    P.run(statement)

###########################
##### merge/aggregate #####
###########################

@follows(mkdir("final_genesets.dir"), assembleWithStringTie)
@merge([assembleWithStringTie,
       os.path.join(
           PARAMS["annotations_dir"],
           PARAMS["annotations_interface_geneset_all_gtf"])],
       "final_genesets.dir/agg-agg-agg.gtf.gz")
def mergeAllAssemblies(infiles, outfile):

    infiles = ["<(zcat %s)" % infile for infile in infiles]
    infiles, reference = infiles[:-1], infiles[-1]

    job_threads = PARAMS["stringtie_merge_threads"]
    job_memory = PARAMS["stringtie_merge_memory"]

    infiles = " ".join(infiles)

    statement = '''stringtie --merge
                            -G %(reference)s
                            -p %(stringtie_merge_threads)s
                            %(stringtie_merge_options)s
                            %(infiles)s
                            2> %(outfile)s.log
                    | cgat gtf2gtf --method=sort
                            --sort-order=gene+transcript
                            -S %(outfile)s -L %(outfile)s.log'''

    P.run(statement) 

###############################################
##### configurable: merge by tissue/group #####
###############################################

@follows(assembleWithStringTie)
@collate([glob.glob("assembled_transcripts.dir/*%s*.gtf.gz" % PARAMS["stringtie_groups"][x]) for x in range(0, len(PARAMS["stringtie_groups"]))],
         regex("(.+)/(.+)-(.+).gtf.gz"),
         add_inputs(os.path.join(
            PARAMS["annotations_dir"],
            PARAMS["annotations_interface_geneset_all_gtf"])),
          r"final_genesets.dir/\2-agg-agg-agg.gtf.gz")
def merge_by_tissue(infiles, outfile):
    job_threads = PARAMS["stringtie_merge_threads"]
    job_memory= PARAMS["stringtie_merge_memory"]

    reference = "<(zcat %s)" % infiles[0][1]
    infiles = ["<(zcat %s)" % infile for infile in infiles[0][0]]
    infiles = " ".join(infiles)
    
    statement = '''stringtie --merge
                    -G %(reference)s
                    -p %(stringtie_merge_threads)s
                    %(stringtie_merge_options)s
                    %(infiles)s
                    2> %(outfile)s.log
                | cgat gtf2gtf --method=sort
                    --sort-order=gene+transcript
                    -S %(outfile)s -L %(outfile)s.log'''
    P.run(statement)

###################
##### utility #####
###################

@follows(mergeAllAssemblies, merge_by_tissue)
def Assembly():
    pass

#######################################################
##### export assembled and aggregated transcripts #####
#######################################################

@follows(mkdir("export/indexed_gtfs.dir"))
@transform([assembleWithStringTie,
            mergeAllAssemblies,
            merge_by_tissue,
            "final_genesets.dir/*.gtf.gz"],
           regex("(.+)/(.+).gtf.gz"),
           r"export/indexed_gtfs.dir/\2.gtf.gz")
def exportIndexedGTFs(infile, outfile):

    statement = ''' zcat %(infile)s
                    | sort -k1,1 -k4,4n
                    | bgzip > %(outfile)s &&
                    tabix -f -p gff %(outfile)s &&
                    if [[ %(outfile)s == *"agg-agg-agg"* ]]; then 
                        cp %(outfile)s %(outfile)s.tbi export/; 
                    fi;
                    '''
    P.run(statement)

###################
##### utility #####
###################

@follows(exportIndexedGTFs)
def export():
    pass

##################################################
##### create decoy-aware transcriptome/index #####
##################################################

@follows(mkdir("salmon_index"), mkdir("salmon_index/DAT"),exportIndexedGTFs)
@transform("final_genesets.dir/agg-agg-agg.gtf.gz",
           formatter(),
           "salmon_index/{basename[0]}.salmon.index")
def makeSalmonIndex(infile,outfile):
    '''Create decoy-aware transcriptome index'''
    job_memory="64G"
    job_threads=1

    gtf_basename = P.snip(os.path.basename(infile), ".gtf.gz")
    transcript_fasta = "salmon_index/" + gtf_basename + "transcripts.fa"
    fastaref =PARAMS["portcullis_fastaref"]
    index_options=PARAMS["salmon_indexoptions"]
    tmpfile = P.get_temp_filename()

    statement = ''' gunzip -c %(infile)s > %(tmpfile)s;
                    gffread %(tmpfile)s -g %(fastaref)s -w %(transcript_fasta)s;
                    grep "^>" <%(fastaref)s | cut -d " " -f 1 > salmon_index/DAT/decoys.txt;
                    sed -i.bak -e 's/>//g' salmon_index/DAT/decoys.txt; 
                    cat %(transcript_fasta)s %(fastaref)s > salmon_index/DAT/gentrome.fa.gz;
                    salmon index
                        -p %(job_threads)s
                        %(index_options)s
                        -t salmon_index/DAT/gentrome.fa.gz
                        -d salmon_index/DAT/decoys.txt
                        -i %(outfile)s;
                    rm %(tmpfile)s
                    '''
    P.run(statement)

##########################
##### quantification #####
##########################

if not os.path.exists("temp_bams/"):
    os.makedirs("temp_bams/")

@follows(mkdir("quantification.dir"), mkdir("sorted_bams"),
         mergeAllAssemblies, makeSalmonIndex)
@product(["input_assemble.dir/*.bam",
          "input_assemble.dir/*.remote"],
         formatter(".+/(?P<TRACK>.+).([bam|remote])"),
         mergeAllAssemblies,
         formatter(".+/(?P<GENESET>.+).gtf.gz"),
         "quantification.dir/{TRACK[0][0]}_{GENESET[1][0]}")
def quantifyWithSalmon(infiles, outfile):
    '''Quantify existing samples against genesets'''
    job_threads=2
    job_memory="24G"

    infile, gtffile = infiles
    basefile = os.path.basename(infile)
    sample_name = basefile.split(os.extsep, 1)
    gtfbase = P.snip(os.path.basename(gtffile), ".gz")
    salmonIndex = "salmon_index/" + gtfbase + ".salmon.index"

    salmon_options = PARAMS["salmon_quantoptions"]
    bootstrap_options = PARAMS["num_bootstraps"]

    sorted_bam="sorted_bams/" + sample_name[0] + "_sorted.bam"
    fastq1 = P.snip(outfile, "_agg-agg-agg")+".1.fastq"
    fastq2 = P.snip(outfile, "_agg-agg-agg")+".2.fastq"
    fastq0 = P.snip(outfile, "_agg-agg-agg")+".0.fastq"

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
                            -o %(outfile)s
                            %(salmon_options)s
                            --numBootstraps %(bootstrap_options)s;
                    else
                        salmon quant -i %(salmonIndex)s
                            --libType IU
                            -1 %(fastq1)s
                            -2 %(fastq2)s
                            -o %(outfile)s
                            %(salmon_options)s
                            --numBootstraps %(bootstrap_options)s;
                    fi; 
                    mv %(outfile)s/quant.sf %(outfile)s.sf; 
                    rm %(fastq1)s; rm %(fastq2)s; rm %(fastq0)s; rm %(sorted_bam)s 
                    '''
    P.run(statement)

#############################################
##### load quantification into database #####
#############################################

@follows(quantifyWithSalmon)
@merge("quantification.dir/*.sf", "database_load/salmon_quant.load")
def mergeAllQuants(infiles, outfile):
    job_threads=3

    P.concatenate_and_load(infiles, outfile,
                         regex_filename="quantification.dir/(.*)_agg-agg-agg.sf",
                         options="-i Name -i Length -i EffectiveLength"
                         " -i TPM -i NumReads -i track"
                         " -i source",
                         job_memory="64G")

    if not os.path.isfile("mapping_rates.txt"):
        statement = ''' bash /shared/sudlab1/General/projects/UTRONs/MyFiles/scripts/mapping_rates_script.sh '''
        P.run(statement)
    else:
        pass

###################################
##### Export data to database #####
###################################

@follows(mergeAllQuants, mkdir("expression.dir", "expression.dir/csvdb_files"))
def CSVDBfiles():
    ''' utility function to connect to database. Use this method to connect to the pipeline database.
        Export all_utrons_ids, novel_utrons_ids and tx2gene data in .txt format to be used in the Rscript. 
        '''

    subprocess.call(["sqlite3", PARAMS["database_name"],
                     ".headers on", ".mode tab", ".output expression.dir/csvdb_files/tx2gene.txt",
                     "select transcript_id, match_gene_id from transcript_class where track = 'agg-agg-agg'"])

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


#def build3UTRlib():


#def extract3UTRseq():


#def quant3UTRusage():


#def compareQAPA():


##################
###### misc ######
##################

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))