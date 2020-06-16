#!/bin/bash

# How many cores can we use
CORES=8

# Prepare a fresh working directory
DIRECTORY=testFolder
if [ -d $DIRECTORY ]; then
  rm -rf $DIRECTORY
fi

./bin/simulator ./data/lambda_virus.fa 150 100 0.01 lambda_virus.reads;

# Start BWise (it creates $DIRECTORY)
./Bwise  -u lambda_virus.reads.fa  -t $CORES  -k 31 -K 31 -o testFolder

# Test ok?
if [ -f "$DIRECTORY/contigs_k31.fa" ];
then
  echo "IT WORKS !";
  ./src/n50 $DIRECTORY/contigs_k31.fa;
else
   echo "FAIL"
fi
