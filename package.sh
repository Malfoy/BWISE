#!/bin/bash



#THIS INSTALLATION REQUIRES GCC, GIT, MAKE, CMAKE



function help {
echo "BWISE installation script"
echo "This installation requires GCC>=4.9, GIT, MAKE and CMAKE3 and Python3"
echo "-f absolute path of folder to put the binaries"
echo "-t to use multiple thread for compilation (default 8)"
}





threadNumber=2
#SCRIPT=$(readlink -f $0)

# Absolute path this script is in. /home/user/bin
folder=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#folder=`dirname $SCRIPT`
folder+="/bin"



cd src;

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
t)
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



mkdir -p $folder 2>/dev/null;
if [ $? -ne 0 ]
       then
              echo "there was a problem with "$folder$" directory creation"
              exit 1
       fi
echo "I put binaries in $folder";



make -j $threadNumber;
if [ $? -ne 0 ]
       then
              echo "there was a problem with binary compilation, check logs"
              exit 1
       fi
# cp bwise ..;
cp K2000/*.py $folder;
cp K2000/run_K2000.sh $folder;
cp sequencesToNumbers $folder;
cp numbersFilter $folder;
cp numbersToSequences $folder;
echo PHASE ZERO LAUNCHER: BWISE;


git clone --recursive https://github.com/GATB/bloocoo.git ;
cd bloocoo;
mkdir build32 2>/dev/null; cd build32;
cmake -DKSIZE_LIST="32" .. ;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bloocoo32 cmake, check logs"
              exit 1
       fi
make -j $threadNumber ;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bloocoo32 compilation, check logs"
              exit 1
       fi
cp bin/Bloocoo Bloocoo32;
cp Bloocoo32 $folder;
cd ..;
mkdir build64 2>/dev/null; cd build64;
cmake -DKSIZE_LIST="64" ..;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bloocoo64 cmake, check logs"
              exit 1
       fi
make -j $threadNumber ;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bloocoo64 compilation, check logs"
              exit 1
       fi
cp bin/Bloocoo Bloocoo64;
cp Bloocoo64 $folder;
cd ../..;
#~ mkdir build128 2>/dev/null; cd build128;
#~ cmake -DKSIZE_LIST="128" .. ;
#~ if [ $? -ne 0 ]
       #~ then
              #~ echo "there was a problem with bloocoo128 cmake, check logs"
              #~ exit 1
       #~ fi
#~ make -j $threadNumber;
#~ if [ $? -ne 0 ]
       #~ then
              #~ echo "there was a problem with bloocoo128 compilation, check logs"
              #~ exit 1
       #~ fi
#~ cp bin/Bloocoo Bloocoo128;
#~ cp Bloocoo128 $folder;
#~ cd ../..;
#~ cp bloocoo/build32/ext/gatb-core/bin/h5dump $folder;
echo PHASE ONE, READ CORRECTION: BLOOCOO;



git clone --recursive https://github.com/GATB/bcalm ;
cd bcalm;
mkdir build 2>/dev/null; cd build;
cmake -DKSIZE_LIST="32 64 128 320 416 1024"  ..  ;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bcalm cmake, check logs"
              exit 1
       fi
make -j $threadNumber ;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bcalm compilation, check logs"
              exit 1
       fi
cp bcalm $folder;
cd ../..;

echo PHASE TWO, GRAPH CONSTRUCTION: BCALM;




git clone https://github.com/Malfoy/BGREAT2 ;
cd BGREAT2;
make -j $threadNumber >>logCompile;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bgreat compilation, check logs"
              exit 1
       fi
cp bgreat $folder;
cd ..;
echo PHASE THREE, READ MAPPING ON THE DBG: BGREAT;




git clone https://github.com/Malfoy/BTRIM ;
cd BTRIM;
make -j $threadNumber;
if [ $? -ne 0 ]
       then
              echo "there was a problem with btrim compilation, check logs"
              exit 1
       fi
cp btrim $folder;
cd ..;
echo PHASE FOUR GRAPH CLEANING: BTRIM;



cd .. ;
tar -czvf bin.tar.gz bin ;



echo The end !;

