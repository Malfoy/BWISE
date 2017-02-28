BWISE
=====

de Bruijn Workflow using Integral information of Short paired-End reads

Work in progress -- so far, this assembler was only tested with > 100X of 2*250 bp Illumina reads (recommended fragment library size ~500-600 bp). 

INSTALLATION
------------

`./install2.sh -f /home/malfoy/BWISEbin` will download and compile all the needeed programs (requires GCC\>=4.9, GIT, MAKE and CMAKE3).

-f absolute path of folder to put the binaries

-t use multiple threads for compilation (default is 8)

You can test your installation by running `./test.sh`

RUN
---

`./bwise -x examplePairedReads.fa -o workingFolder`

Options and default values:

-x paired read file (at most one at present, with the reads interleaved)

-u unpaired read file (at most one at present)

-o working folder (Bwise will run the complete pipeline in this folder, default is .)

-s kmer solidity threshold (only kmers occurring at least s times are included in the de Bruijn graph; default is 2)

-S unitig solidity threshold (2)

-k maximal kmer size (for de Bruijn graph construction; default is 240)

-p superReads cleaning threshold (2)

-c number of correction steps (default is at most 3)

-t maximum number of threads 

 
