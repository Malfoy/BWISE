#!/bin/bash



#THIS INSTALLATION REQUIRES GCC, GIT, MAKE, CMAKE



function help {
echo "BWISE installation script"
echo "This installation requires GCC>=4.9, GIT, MAKE and CMAKE3 and Python3"
echo "-p absolute path of folder to put the binaries"
echo "-t to use multiple thread for compilation (default 8)"
echo "-f Fast install, download the bcalm binaries directly, do not need CMAKE3 but large k (>101) are unavailable"
}



threadNumber=8
FASTINSTALL=0

# Absolute path this script is in. /home/user/bin
folderInstall=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
folder=$folderInstall"/bin"



while getopts "hfp:t:" opt; do
case $opt in
f)
echo "FAST  installation:"
FASTINSTALL=1
;;
h)
help
exit
;;
p)
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




echo PHASE ONE LAUNCHER: BWISE;


cat src/Bwise_header.py>Bwise
echo "BWISE_MAIN = os.path.dirname(\"$folder\")" >> Bwise
cat src/Bwise_broken.py>>Bwise
cd src;


make LOL=-Dfolder=$folder -j $threadNumber >>logCompile.txt 2>>logCompile.txt;
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
cp path_to_kmer $folder;


echo PHASE TWO, GRAPH CONSTRUCTION: BCALM;


if [ $FASTINSTALL -eq 1 ]; then
	echo "fast install of bcalm"
	wget https://github.com/GATB/bcalm/releases/download/v2.2.3/bcalm-binaries-v2.2.3-Linux.tar.gz >>logCompile.txt 2>>logCompile.txt;

	tar -xvf bcalm-binaries-v2.2.3-Linux.tar.gz >>logCompile.txt 2>>logCompile.txt;

	cp bcalm-binaries-v2.2.3-Linux/bin/bcalm $folder;
fi


if [ $FASTINSTALL -eq 0 ]; then
	echo "bcalm compilation"
git clone --recursive https://github.com/GATB/bcalm --depth 1 >>logCompile.txt 2>>logCompile.txt;
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
fi



echo PHASE THREE, READ MAPPING ON THE DBG: BGREAT;


git clone https://github.com/Malfoy/BGREAT2 --depth 1 >>logCompile.txt 2>>logCompile.txt;
# git checkout 6a5afe388ccf733a3c73ff3f9d912174a0697fa8>>../../logCompile.txt 2>>../../logCompile.txt;

cd BGREAT3;
make -j $threadNumber >>../logCompile.txt 2>>../logCompile.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with bgreat compilation, check logs"
              exit 1
       fi
cp bgreat $folder;
cd ..;




echo PHASE FOUR GRAPH CLEANING: BTRIM;
	


git clone https://github.com/Malfoy/BTRIM --depth 1 >>logCompile.txt 2>>logCompile.txt;


cd BTRIM;
make -j $threadNumber >>../logCompile.txt 2>>../logCompile.txt.txt;
if [ $? -ne 0 ]
       then
              echo "there was a problem with btrim compilation, check logs"
              exit 1
       fi
cp btrim $folder;
cd ..;





echo The end !;
