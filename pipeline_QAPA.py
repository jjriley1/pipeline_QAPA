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



Input files
===========



Requirements
============



Pipeline output
===============



Code
====

"""