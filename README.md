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

-h, --help            show this help message and exit
  -x PAIRED_READFILES   input fasta or (compressed .gz if -c option is != 0)
                        paired-end read files. Several read files must be
                        concatenated.
  -u SINGLE_READFILES   input fasta or (compressed .gz if -c option is != 0)
                        single-end read files. Several read files must be
                        concatenated.
  -c NB_CORRECTION      an integer, number of steps of read correction
                        (default 1)
  -s KMER_SOLIDITY      an integer, k-mers present strictly less than this
                        number of times in the dataset will be discarded
                        (default 2)
  -S KMER_COVERAGE      an integer, minimal unitig coverage for first cleaning
                        (default 5)
  -p SR_SOLIDITY        an integer, super-reads present strictly less than
                        this number of times will be discarded (default 2)
  -P SR_COVERAGE        an integer X, unitigs with less than size/X reads
                        mapped is filtred (default 20)
  -k K_MIN              an integer, smallest k-mer size (default 63)
  -K K_MAX              an integer, largest k-mer size (default 201)
  -e MAPPING_EFFORT     Anchors to test for mapping (default 100)
  -a ANCHOR_SIZE        Anchors size (default 41)
  -m MISSMATCH_ALLOWED  missmatch allowed in mapping (default 2)
  -t NB_CORES           number of cores used (default max)
  -o OUT_DIR            path to store the results (default = current
                        directory)
  --version             show program's version number and exit
