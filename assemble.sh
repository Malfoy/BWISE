#!/bin/bash



function help {
echo "BWISE"
echo "-x for paired read file"
echo "-u for unpaired read file"
echo "-k for kmer size"
echo "-o for working folder"
echo "-x for paired read file"
echo "-s for kmer solidity threshols"
}



#ARGUMENT PARSING


outputFolder="folderAssembly"
readFileName=""
pairedReadFileName=""
kmersize=100
solidity=5

while getopts "hx:u:k:o:s:" opt; do
case $opt in

h)
help
exit
;;

x)
echo "use paired read file: $OPTARG" >&2
pairedReadFileName=$OPTARG
;;

s)
echo "use kmer solidity threshold: $OPTARG" >&2
solidity=$OPTARG
;;

u)echo "use unpaired read file: $OPTARG" >&2
readFileName=$OPTARG
;;

o)echo "use output folder: $OPTARG" >&2
outputFolder=$OPTARG
;;


k)
echo "use k: $OPTARG" >&2
kmersize=$OPTARG
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



#CLEAN DIR



mkdir $outputFolder;
cd $outputFolder;
rm *>>log 2>>log;
start_time=`date +%s`



if [ -z "$readFileName"  ]; then
	if [ -z "$pairedReadFileName"  ]; then
	help
	exit 0
	fi
fi

bloocooString="../$pairedReadFileName,../$readFileName"
echo "reads_corrected_0_.fasta
reads_corrected_1_.fasta">bankBcalm;
bgreatString="-x reads_corrected_0_.fasta -u reads_corrected_1_.fasta"

if [ -z "$readFileName" ]; then
	echo "case 1"
	bloocooString="../$pairedReadFileName"
	rm bankBcalm;
	echo "reads_corrected.fa" > bankBcalm;
	bgreatString="-x reads_corrected.fa"
fi

if [ -z "$pairedReadFileName" ]; then
	echo "case 2"
	bloocooString="../$readFileName"
        rm bankBcalm;
	echo "reads_corrected.fa" > bankBcalm;
	bgreatString="-u reads_corrected.fa"
fi



#READ CORRECTION



startcorrection_time=`date +%s`
../bloocoo//build/bin/Bloocoo -file $bloocooString  -kmer-size 31 -abundance-min 5 -out reads_corrected.fa >>log 2>>log;
endcorrection_time=`date +%s`
echo 1/5 reads corrected in `expr $endcorrection_time - $startcorrection_time` seconds;



#DBG CONSTRUCTION



startgraph_time=`date +%s`
../bcalm/build/bcalm -in bankBcalm -kmer-size $kmersize -abundance-min $solidity -out out >>log 2>>log;
../kMILL/src/kMILL out.unitigs.fa $((kmersize-1)) $((kmersize-2))>>log 2>>log ;
mv out_k$((kmersize-1))_.fa out.fa;
endgraph_time=`date +%s`
echo 2/5 graph constructed in `expr $endgraph_time - $startgraph_time` seconds;



#READ MAPPING ON DBG



startmap_time=`date +%s`
../BGREAT2/bgreat -k $kmersize $bgreatString -g out.fa -t 8  -c -m 0 -e 1 >>log 2>>log;
endmap_time=`date +%s`
echo 3/5 read mapped on graph in `expr $endmap_time - $startmap_time` seconds;



#SUPERREADS CLEANING



#duplicate superReads elimination
startclean_time=`date +%s`
../kMILL/src/pathsCleaner paths 3	 >>log 2>>log;
#elimination of contained superreads (non maximal one)
echo "noduplicate.fa" > bank
../BREADY/short_read_connector.sh -b bank -q bank >> log 2>>log;
rm *.h5;
endclean_time=`date +%s`
echo 4/5 paths redundancy cleaned in `expr $endclean_time - $startclean_time` seconds;



#COMPACTION OF SUPERREADS



startass_time=`date +%s`
../kMILL/src/kMILL short_read_connector_res.txt 50000 >>log 2>>log;
endass_time=`date +%s`
echo 5/5 maximal paths compacted in `expr $endass_time - $startass_time` seconds;
mv out_k50000_.fa contigs1.fa;



#END



end_time=`date +%s`
echo Total execution time was `expr $end_time - $start_time` seconds;
#number of contigs obtained
raw_n50 contigs1.fa;




#~ rm noduplicate.fa notAligned.fa paths reads_corrected.fasta reads_corrected.unitigs.fa reads.fasta short_read_connector_res.txt >>log 2>>log;



