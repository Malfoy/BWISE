BWISE
=====

de Bruijn Workflow using Integral information of Short paired-End reads

Work in progress -- so far, this assembler was only tested with > 100X of 2*250 bp Illumina reads (recommended fragment library size ~500-600 bp). 

INSTALLATION
------------

`./install2.sh -f /home/malfoy/BWISEbin`

Will download and compile all the needeed programs.

 

This installation requires GCC\>=4.9, GIT, MAKE and CMAKE3
-f absolute path of folder to put the binaries
-t to use multiple threads for compilation (default is 8)

You can test your installation  with `./test.sh`

RUN
---

`./bwise -x examplePairedReads.fa -o workingFolder`

Bwise will run the complete pipeline in the directory folderAssembly 
-x for paired read file 
-u for unpaired read file 
-o for working folder 
-s for kmer solidity threshold 
-k for largest kmer size 
-p for superReads cleaning threshold
