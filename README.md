BWISE
=====

de Bruijn Workflow using Integral information of Short paired-End reads

Work in progress -- so far, this assembler works well  with \> 100X of
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


UPDATE
------------

`./update.sh` will update the components of Bwise

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

  -x PAIRED_READFILE    input fasta/fastq (fa, fa.gz, fq or fq.gz) paired-end read files. Only one input file of this type can be used: if you have for instance two paired-end librairies (e.g. libA and libB), first interleave the forward and reverse reads for each library (fq2fa --merge libA_1.fq libA_2.fq interleavedlibA.fa; fq2fa --merge libB_1.f libB_2.fq interleavedlibB.fa using for instance the tool fq2fa of the IDBA package), then concatenate the two resulting files (cat interleavedlibA.fa interleavedlibB.fa > interleavedlibsA+B.fa) in order to generate the PAIRED_READFILE input fasta).

  -u SINGLE_READFILE    input fasta/fastq (fa, fa.gz, fq or fq.gz) single-read files. Only one input file of this type can be used: if you have several single-read librairies, concatenate them to generate the SINGLE_READFILE input fasta).


  -s KMER_SOLIDITY      an integer, k-mers present strictly less than this
                        number of times in the dataset will be discarded
                        (default 2)

  -S KMER_COVERAGE      an integer, minimal unitig coverage for first cleaning
                        (default auto)

  -p SR_SOLIDITY        an integer, super-reads present strictly less than
                        this number of times will be discarded (default 2)

  -P SR_COVERAGE        an integer X, unitigs with less than X reads
                        mapped is filtred (default 10)

  -k K_MIN              an integer, smallest k-mer size (default 63)

  -K K_MAX              an integer, largest k-mer size (default 201). Advice: put a -K near slightly below your read size (241 for 250 bp reads for example)

  -e MAPPING_EFFORT     number of anchors to test for mapping (default max)

  -a ANCHOR_SIZE        size of the anchors (default 31)

  -m MISSMATCH_ALLOWED  number of missmatchs allowed in mapping (default 10)

  -t NB_CORES           number of cores used (default 1)

  -o OUT_DIR            path to store the results (default = current
                        directory)

  --version             show program's version number and exit
