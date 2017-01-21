#!/bin/bash

mkdir folderTest;

./bwise -x examplePairedReads.fa -u exampleUnpairedReads.fa  -o folderTest -c 4;

if [ -f "folderTest/contigs.fa" ];
then
	echo "IT WORKS !";
	./n50 folderTest/contigs.fa;
else
   echo "FAIL"
fi


