#!/bin/bash



#THIS INSTALLATION REQUIRES GCC, GIT, MAKE, CMAKE



function help {
echo "BWISE installation script"
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



cat src/Bwise_header.py>Bwise.py
echo "BWISE_MAIN = os.path.dirname(\"$folder\")" >> Bwise.py
cat src/Bwise_broken.py>>Bwise.py
cd src;


make LOL=-Dfolder=$folder -j $threadNumber >>logCompile.txt 2>>logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with binary compilation, check logs"
              exit 1
       fi
# cp bwise ..;
cp K2000/*.py $folder;
cp K2000/*.sh $folder;
cp sequencesToNumbers $folder;
cp numbersFilter $folder;
cp numbersToSequences $folder;
cp simulator $folder;



echo PHASE ONE LAUNCHER: BWISE;



git clone --recursive https://github.com/GATB/bcalm >>logCompile.txt 2>>logCompile.txt;
cd bcalm;
mkdir build 2>/dev/null; cd build;
cmake -DKSIZE_LIST="32 64 128 256 512 1024"  ..  >>../logCompile.txt 2>>../logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bcalm cmake, check logs"
              exit 1
       fi
make -j $threadNumber >>../logCompile.txt 2>>../logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bcalm compilation, check logs"
              exit 1
       fi
cp bcalm $folder;
cd ../..;



echo PHASE TWO, GRAPH CONSTRUCTION: BCALM;



git clone https://github.com/Malfoy/BGREAT2 >>logCompile.txt 2>>logCompile.txt;
cd BGREAT2;
make -j $threadNumber >>../logCompile.txt 2>>../logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bgreat compilation, check logs"
              exit 1
       fi
cp bgreat $folder;
cd ..;



echo PHASE THREE, READ MAPPING ON THE DBG: BGREAT;



git clone https://github.com/Malfoy/BTRIM >>logCompile.txt 2>>logCompile.txt;
cd BTRIM;
make -j $threadNumber >>../logCompile.txt 2>>../logCompile.txt.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with btrim compilation, check logs"
              exit 1
       fi
cp btrim $folder;
cd ..;



echo PHASE FOUR GRAPH CLEANING: BTRIM;



echo The end !;

