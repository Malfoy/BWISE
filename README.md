BWISE 

de Bruijn Workflow using Integral information of Short pair End reads

INSTALATION

./install2.sh -f /home/malfoy/BWISEbin
 
Will download and compile all the needeed software
This instalation need GCC>=4.9, GIT, MAKE, CMAKE3 !
-f absolute path of folder to put the binaries 
-t to use multiple thread for compilation (default 8)

You can test your install with 
./test.sh


RUN

./bwise -x examplePairedReads.fa -o workingFolder

Bwise will run the complete pipeline in the directory folderAssembly
-x for paired read file
-u for unpaired read file
-o for working folder
-s for kmer solidity threshold
-k for largest kmer size
-p for superReads cleaning threshold

