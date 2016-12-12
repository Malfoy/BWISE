#!/bin/bash



#THIS INSTALATION NEED GCC, GIT, MAKE, CMAKE



function help {
echo "BWISE instalation script"
echo "This instalation need GCC>=4.9, GIT, MAKE, CMAKE3"
echo "-f absolute path of folder to put the binaries"
echo "-t to use multiple thread for compilation (default 8)"
}



threadNumber=8
folder=""


while getopts "hf:t:" opt; do
case $opt in
h)
help
exit
;;
f)
echo "use folder: $OPTARG" >&2
folder=$OPTARG
;;
s)
echo "use  $OPTARG threads" >&2
threadNumber=$OPTARG
;;
\?)
echo "Invalid option: -$OPTARG" >&2
exit 1
;;
:)
echo "Option -$OPTARG requires an argument." >&2
exit 1
;;
esac
done

if [ -z "$folder"  ]; then
	help
	exit 0
	fi



mkdir $folder;



echo PHASE ZERO LAUNCHER: BWISE;
make LOL=-Dfolder=$folder -j $theadNumber >>logCompile 2>>logCompile;
cp bwise $folder;



echo PHASE ONE, READ CORRECTION: BLOOCOO;
git clone --recursive https://github.com/GATB/bloocoo.git >>logCompile 2>>logCompile;
cd bloocoo;
mkdir build; cd build;
#TODO multiple bin
cmake  .. >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bin/Bloocoo $folder;
cd ../..;




echo PHASE TWO, GRAPH CONSTRUCTION: BCALM;
git clone --recursive https://github.com/GATB/bcalm >>logCompile 2>>logCompile;
cd bcalm;
mkdir build; cd build;
cmake .. -DKSIZE_LIST="32 64 96 128 160 192 224 256" >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bcalm $folder;
cd ../..;




echo PHASE THREE, READ MAPPING ON THE DBG: BGREAT;
git clone https://github.com/Malfoy/BGREAT2 >>logCompile 2>>logCompile;
cd BGREAT2;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bgreat $folder;
cd ..;




echo PHASE FOUR, SUPERREADS CLEANING: BREADY;
git clone --recursive https://github.com/Malfoy/BREADY >>logCompile 2>>logCompile;
cd BREADY;
mkdir build; cd build;
cmake .. >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bin/BREADY $folder;
cd ../..;




echo and DSK;
git clone --recursive https://github.com/GATB/dsk.git >>logCompile 2>>logCompile;
cd dsk;
mkdir build;  cd build;
cmake -DKSIZE_LIST="64" .. >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bin/dsk $folder;
cd ../..;



echo PHASE FIVE, MAXIMAL SUPERREADS COMPACTION: KMILL;
git clone https://github.com/kamimrcht/kMILL >>logCompile 2>>logCompile;
cd kMILL/src;
make -j $threadNumber >>logCompile 2>>logCompile;
cp kMILL $folder;
cp pathsCleaner $folder;
cd ../..;


echo The end !;
