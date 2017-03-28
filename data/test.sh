#!/bin/bash

# How many cores can we use
CORES=2

# Prepare a fresh working directory
DIRECTORY=folderTest
if [ -d $DIRECTORY ]; then
  rm -rf $DIRECTORY
fi

# Start BWise (it creates $DIRECTORY)
../bwise -x examplePairedReads.fa -u exampleUnpairedReads.fa  -o $DIRECTORY -c $CORES

# Test ok?
if [ -f "$DIRECTORY/contigs.fa" ];
then
  echo "IT WORKS !";
  ../src/n50 $DIRECTORY/contigs.fa;
else
   echo "FAIL"
fi

