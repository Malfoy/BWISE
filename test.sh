#!/bin/bash

mkdir folderTest;
./bwise -x examplePairedReads.fa -o folderTest >>logtest 2>>logtest ;
if [ -f "folderTest/contigs.fa" ];
then
   echo "IT WORKS !"
else
   echo "FAIL"
fi         


