BWISE
=====

de Bruijn Workflow using Integral information of Short paired-End reads

Work in progress -- so far, this assembler was only tested with \> 100X of
2\*250 bp Illumina reads (recommended fragment library size \~500-600 bp).


[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

[![Build Status](https://travis-ci.org/Malfoy/BWISE.svg?branch=master)](https://travis-ci.org/Malfoy/BWISE)


INSTALLATION
------------

`./install.sh` will download and compile all the needeed programs (requires
GCC\>=4.9, GIT, MAKE and CMAKE3).

Possible compilation options:

\-f absolute path of folder to put the binaries

\-t use multiple threads for compilation (default is 8)

You can test your installation as follows:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ bash
cd data
./test.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RUN
---

`python3 Bwise.py -x examplePairedReads.fa -o workingFolder`

Options and default values:

 \-h, --help            shows help message and exit
 
 \-x paired reads file (at most one file at present, with the read pairs interleaved)
 
 \-u unpaired read file (at most one file at present)
 
 \-s kmer solidity threshold (only kmers occurring at least s times are included in the initial de Bruijn graph; default is 2, which seems a good value when starting from \~100X data, but if the read coverage is very high this can be raised to 5 or even 10
 
 \-S unitig solidity threshold (only unitigs on which at least S/size reads map are considered valid and retained; default is 2)
 
 \-o output folder (Bwise will run the complete pipeline in this folder; default is the current directory)
 
 \-k maximal kmer size (for de Bruijn graph construction; default is 201)
 
 \-p superReads cleaning threshold (only super-reads that appear more than p times are retained for the final super-read assembly step; default is 2)
 
 \-c number of correction steps (default is 2)
 
 \-t maximum number of threads
 
 \-e numbers of anchors tested during mapping (default is 100)
 
 \-C unitig coverage for the first cleaning (default is 5)
 
 \-m number of missmatches allowed during mapping (default is 2)
 
 \--version shows the currently installed version of the program then exits
