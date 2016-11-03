#!/bin/bash



#THIS SCRIPT NEED GCC, GIT, MAKE, CMAKE, PYTHON



#PHASE ONE, READ CORRECTION: BLOOCOO
cd bloocoo;
git pull;
cd build; make -j 8;
cd ../..;




#PHASE TWO, GRAPH CONSTRUCTION: BCALM
cd bcalm;
git pull;
cd build; make -j 8;
cd ../..;




#PHASE THREE, READ MAPPING ON THE DBG: BGREAT
cd BGREAT2;
git pull;
make -j 8;
cd ..;



#PHASE FOUR, SUPERREADS CLEANING: BREADY
cd BREADY;
git pull;
cd build; make -j 8;
cd ../..;



#PHASE FIVE, MAXIMAL SUPERREADS COMPACTION: KMILL
cd kMILL/src;
git pull;
make -j 8;
cd ../..;
