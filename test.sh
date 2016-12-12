#!/bin/bash

mkdir folderTest;
./bwise -x examplePairedReads.fa -o folderTest;
if [ -f "folderTest/contigs.fa" ];
then
   echo "OK"
else
   echo "FAIL"
fi         


