# DUMMY FILE, update with utility functions

'''PipelineGenotypeQC.py- Utility functions for genotype quality control
========================================================================

This is a module file to support the pipeline_genotype_QC.py pipeline.

See https://github.com/CGATOxford/CGATPipelines/blob/master/CGATPipelines/PipelineMapping.py
for an example.
Also see https://www.cgat.org/downloads/public/cgatpipelines/documentation/
for the general tools and approach being followed.

:Release: $Id$
:Date: |today|
:Tags: Python

Genotype QC xxx 

It implements:
   * .bed: plink formatted genotype data
   * 

The basic class :class:`xxx`.  The aim of this
class and its derivatives is to build a sequence of command line
statements that can be sent to a single node on a cluster to process
the input data.

-----
The basic usage inside a pipeline task is as such::
    @transform()
    def mapReads(infile, outfile):
        # initialize the Tool
        m = PipelineMapping.Hisat(
             executable=P.substituteParameters(**locals())["hisat_executable"],
             strip_sequence=PARAMS["strip_sequence"])
        # build the command line statement
        statement = m.build((infile,), outfile)
        P.run()
-----

When implementing a tool, avoid specifying algorithmic options as
class variables. Instead use an option string that can be set in
:file:`pipeline.ini`. The only arguments to a tool constructor should
pertain to pipeline integration, such as filenames, index locations,
threading and in general processing options that change the tools
input/output, as these need to be tracked by the pipeline.
The benefit of this approach is to provide complete control to the
user and is likely to work with different versions of a tool, i.e., if
a command line option to a tool changes, just the configuration file
needs to be changed, but the code remains the same.

Requirements:
* samtools >= 1.1
* plink >=1.90

Reference
---------
'''

import os
import glob
import collections
import re
import gzip
import itertools
import CGATPipelines.Pipeline as P
import logging as L
import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def xxx(xxx):
    '''
    xxx
    Arguments
    ---------
    varxxxtrack : 

    Returns
    -------

    '''
DO SOMETHING
