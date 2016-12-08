#!/bin/bash



#THIS SCRIPT NEED GCC, GIT, MAKE, CMAKE, PYTHON



#PHASE ONE, READ CORRECTION: BLOOCOO
git clone --recursive https://github.com/GATB/bloocoo.git;
cd bloocoo;
mkdir build;  cd build;  cmake DKSIZE_LIST="128"  ..;  make -j 8;
cd ../..;




#PHASE TWO, GRAPH CONSTRUCTION: BCALM
git clone --recursive https://github.com/GATB/bcalm ;
cd bcalm;
mkdir build;  cd build; 
 cmake .. -DKSIZE_LIST="32 64 96 128 160 192 224 256 320";
  make -j 8;
cd ../..;



#PHASE THREE, READ MAPPING ON THE DBG: BGREAT
git clone https://github.com/Malfoy/BGREAT2;
cd BGREAT2;
make -j 8;
cd ..;



#PHASE FOUR, SUPERREADS CLEANING: BREADY
git clone --recursive https://github.com/Malfoy/BREADY;
cd BREADY;
mkdir build;  cd build;  cmake ..;  make -j 8;
cd ../..;



rm BREADY/thirdparty/dsk/bin/linux/*;
git clone --recursive https://github.com/GATB/dsk.git;
cd dsk;
mkdir build;  cd build;  cmake ..;  make -j 8;
cd ../..;
cp dsk/build/bin/* BREADY/thirdparty/dsk/bin/linux/


#PHASE FIVE, MAXIMAL SUPERREADS COMPACTION: KMILL
git clone https://github.com/kamimrcht/kMILL;
cd kMILL/src;
make -j 8;
cd ../..;
