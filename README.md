# pipeline_QAPA

## Overview

This pipeline functions to quantify levels of alternative polyadenylation (APA)
using the QAPA package released by Ha, Blencowe and Morris (2018). The output
from QAPA consists of a database containing individual sample APA quantifications.
Custom python and R scripts are used to interogate this database and subsequently
allow for multiple comparisons to be made based on the design.tsv (see below). 

## Usage


To execute this pipeline in an interactive ShARC session or from a personal
device (not recommended) use the following command:

`python <path_to_pipeline_folder>/pipeline_QAPA.py make full -v5`

To execute this pipeline to the ShARC cluster, use the following:

`submit_pipeline <path_to_pipeline_folder>/pipeline_QAPA.py make full -v5`
