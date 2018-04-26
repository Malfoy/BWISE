#!/bin/bash



#THIS INSTALLATION REQUIRES GCC, GIT, MAKE, CMAKE



function help {
echo "BWISE update script"
echo "This installation requires GCC>=4.9, GIT, MAKE and CMAKE3 and Python3"
echo "-f absolute path of folder to put the binaries"
echo "-t to use multiple thread for compilation (default 8)"
}



threadNumber=8

# Absolute path this script is in. /home/user/bin
folderInstall=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
folder=$folderInstall"/bin"



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

#PULL BWISE
git pull >logCompile.txt 2>>logCompile.txt;


#Update BWISE
cat src/Bwise_header.py>Bwise.py;
echo "BWISE_MAIN = os.path.dirname(\"$folder\")" >> Bwise.py;
cat src/Bwise_broken.py>>Bwise.py
cd src;


make LOL=-Dfolder=$folder -j $threadNumber >>../logCompile.txt 2>>../logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with binary compilation, check logs"
              exit 1
       fi
cp K2000/*.py $folder;
cp K2000/*.sh $folder;
cp sequencesToNumbers $folder;
cp numbersFilter $folder;
cp path_counter  $folder;
cp maximal_sr $folder;
cp simulator $folder;
cp crush_bulle $folder;



echo PHASE ONE LAUNCHER: BWISE;


#PULL Bcalm
cd bcalm;
git pull >>../../logCompile.txt 2>>../../logCompile.txt;
cd build;
cmake -DKSIZE_LIST="32 64 128 256 512 1024"  ..  >>../../../logCompile.txt 2>>../../../logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bcalm cmake, check logs"
              exit 1
       fi
make -j $threadNumber >>../../../logCompile.txt 2>>../../../logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bcalm compilation, check logs"
              exit 1
       fi
cp bcalm $folder;
cd ../..;



echo PHASE TWO, GRAPH CONSTRUCTION: BCALM;


cd BGREAT2;
git checkout 6a5afe388ccf733a3c73ff3f9d912174a0697fa8>>../../logCompile.txt 2>>../../logCompile.txt;


make -j $threadNumber >>../../logCompile.txt 2>>../../logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bgreat compilation, check logs"
              exit 1
       fi
cp bgreat $folder;
cd ..;



echo PHASE THREE, READ MAPPING ON THE DBG: BGREAT;


cd BTRIM;
git pull >>../../logCompile.txt 2>>../../logCompile.txt;
make -j $threadNumber >>../../logCompile.txt 2>>../../logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with btrim compilation, check logs"
              exit 1
       fi
cp btrim $folder;
cd ..;



echo PHASE FOUR GRAPH CLEANING: BTRIM;



echo The end !;

