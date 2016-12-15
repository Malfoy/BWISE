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




eval "
make LOL=-Dfolder=$folder -j $threadNumber >>logCompile 2>>logCompile;
cp bwise $folder;
echo PHASE ZERO LAUNCHER: BWISE;
" ;




eval  "
mkdir bloocoo32; cd bloocoo32;
git clone --recursive https://github.com/GATB/bloocoo.git >>logCompile 2>>logCompile;
cd bloocoo;
mkdir build; cd build;
cmake -DKSIZE_LIST="32" .. >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
mv bin/Bloocoo Bloocoo32;
cp Bloocoo32 $folder;
cd ../../..;

mkdir bloocoo64; cd bloocoo64;
git clone --recursive https://github.com/GATB/bloocoo.git >>logCompile 2>>logCompile;
cd bloocoo;
mkdir build; cd build;
cmake -DKSIZE_LIST="64" .. >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
mv bin/Bloocoo Bloocoo64;
cp Bloocoo64 $folder;
cd ../../..;

mkdir bloocoo128; cd bloocoo128;
git clone --recursive https://github.com/GATB/bloocoo.git >>logCompile 2>>logCompile;
cd bloocoo;
mkdir build; cd build;
cmake -DKSIZE_LIST="128" .. >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
mv bin/Bloocoo Bloocoo128;
cp Bloocoo128 $folder;
cd ../../..;

echo PHASE ONE, READ CORRECTION: BLOOCOO;
" ;




eval  "
git clone --recursive https://github.com/GATB/bcalm >>logCompile 2>>logCompile;
cd bcalm;
mkdir build; cd build;
cmake .. -DKSIZE_LIST="32 64 96 128 160 192 224 256" >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bcalm $folder;
cd ../..;
echo PHASE TWO, GRAPH CONSTRUCTION: BCALM;
" ;




eval  "
git clone https://github.com/Malfoy/BGREAT2 >>logCompile 2>>logCompile;
cd BGREAT2;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bgreat $folder;
cd ..;
echo PHASE THREE, READ MAPPING ON THE DBG: BGREAT;
" ;



eval  "

git clone --recursive https://github.com/Malfoy/BREADY >>logCompile 2>>logCompile;
cd BREADY;
mkdir build; cd build;
cmake .. >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bin/BREADY $folder;
cd ../..;
git clone --recursive https://github.com/GATB/dsk.git >>logCompile 2>>logCompile;
cd dsk;
mkdir build;  cd build;
cmake -DKSIZE_LIST="32" .. >>logCompile 2>>logCompile;
make -j $threadNumber >>logCompile 2>>logCompile;
cp bin/dsk $folder;
cd ../..;
echo PHASE FOUR, SUPERREADS CLEANING: BREADY;
" ;



eval  "

git clone https://github.com/kamimrcht/kMILL >>logCompile 2>>logCompile;
cd kMILL/src;
make -j $threadNumber >>logCompile 2>>logCompile;
cp kMILL $folder;
cp sequencesCleaner $folder;
cd ../..;
echo PHASE FIVE, MAXIMAL SUPERREADS COMPACTION: KMILL;

" ;

wait;
echo The end !;
