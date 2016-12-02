#!/bin/bash



function help {
echo "BWISE"
echo "-x for paired read file"
echo "-u for unpaired read file"
echo "-k for kmer size"
echo "-o for working folder"
echo "-s for kmer solidity threshold"
echo "-K for large kmer size"
echo "-S for large solidity threshold"
echo "-t for thread number"
}



#ARGUMENT PARSING


outputFolder="folderAssembly"
readFileName=""
pairedReadFileName=""
kmersize=100
kmersize2=200
solidity=5
solidity2=2
pathsSolidity=2

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

S)
echo "use large kmer solidity threshold: $OPTARG" >&2
solidity2=$OPTARG
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

K)
echo "use K: $OPTARG" >&2
kmersize2=$OPTARG
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



echo "Parameters used for this assembly :
unpaired read file used $readFileName=
paired file used: $pairedReadFileName
kmersize small graph: $kmersize
kmer size large graph: $kmersize2
solidity small graph: $solidity
solidity large graph: $solidity2
path solidity: $pathsSolidity
" > ParametersUsed.txt 


#READ CORRECTION



qqqstartcorrection_time=`date +%s`
../bloocoo/build/bin/Bloocoo -file $bloocooString  -kmer-size 31 -abundance-min 5 -out reads_corrected.fa >>log 2>>log;
# $bloocooString reads_corrected.fa 
endcorrection_time=`date +%s`
echo 1/5 reads corrected in `expr $endcorrection_time - $startcorrection_time` seconds;



#DBG CONSTRUCTION



startgraph_time=`date +%s`
../bcalm/build/bcalm -in bankBcalm -kmer-size $kmersize -abundance-min $solidity -out out >>log 2>>log;
../kMILL/src/kMILL out.unitigs.fa $((kmersize-1)) $((kmersize-2))>>log 2>>log ;
mv out_out.unitigs.fa.fa out.fa;
endgraph_time=`date +%s`
echo 2/5 graph constructed in `expr $endgraph_time - $startgraph_time` seconds;



#READ MAPPING ON DBG



startmap_time=`date +%s`
../BGREAT2/bgreat -k $kmersize $bgreatString -g out.fa -t 20  -c -m 0 -e 1 ;
endmap_time=`date +%s`
echo 3/5 read mapped on graph in `expr $endmap_time - $startmap_time` seconds;



#DBG CONSTRUCTION2



startgraph_time=`date +%s`
../bcalm/build/bcalm -in paths -kmer-size $kmersize2 -abundance-min $solidity2 -out out2 >>log 2>>log;
../kMILL/src/kMILL out2.unitigs.fa $((kmersize2-1)) $((kmersize2-2))>>log 2>>log ;
mv out_out2.unitigs.fa.fa out2.fa;
endgraph_time=`date +%s`
echo 4/5 large graph constructed in `expr $endgraph_time - $startgraph_time` seconds;



#READ MAPPING ON DBG2



startmap_time=`date +%s`
rm paths;
../BGREAT2/bgreat -k $kmersize2 $bgreatString -g out2.fa -t 20  -c -m 0 -e 1 ;
endmap_time=`date +%s`
echo 4/5 read mapped on large  graph in `expr $endmap_time - $startmap_time` seconds;



#SUPERREADS CLEANING



#duplicate superReads elimination
startclean_time=`date +%s`
../kMILL/src/pathsCleaner paths $pathsSolidity	 >>log 2>>log;
#elimination of contained superreads (non maximal one)
echo "noduplicate.fa" > bank
../BREADY/short_read_connector.sh -b bank -q bank >> log 2>>log;
rm *.h5;
endclean_time=`date +%s`
echo 5/5 paths redundancy cleaned in `expr $endclean_time - $startclean_time` seconds;



#COMPACTION OF SUPERREADS



startass_time=`date +%s`
../kMILL/src/kMILL short_read_connector_res.txt >>log 2>>log;
endass_time=`date +%s`
echo 6/5 maximal paths compacted in `expr $endass_time - $startass_time` seconds;
mv out_short_read_connector_res.txt.fa contigs1.fa;



#END



end_time=`date +%s`
echo Total execution time was `expr $end_time - $start_time` seconds;
#number of contigs obtained
n50 contigs1.fa;



exit 0;






startgraph_time=`date +%s`
../bcalm/build/bcalm -in contigs1.fa -kmer-size 200 -abundance-min 1 -out out2 >>log 2>>log;
../kMILL/src/kMILL out2.unitigs.fa 199 198>>log 2>>log ;
mv out_out2.unitigs.fa.fa out2.fa;
endgraph_time=`date +%s`
echo 6/5 graph constructed in `expr $endgraph_time - $startgraph_time` seconds;



startmap_time=`date +%s`
../BGREAT2/bgreat -k 200 $bgreatString -g out2.fa -t 20  -c -m 0 -e 1 ;
endmap_time=`date +%s`



#duplicate superReads elimination
startclean_time=`date +%s`
../kMILL/src/pathsCleaner paths 3   >>log 2>>log;
#elimination of contained superreads (non maximal one)
echo "noduplicate.fa" > bank
../BREADY/short_read_connector.sh -b bank -q bank >> log 2>>log;
rm *.h5;
endclean_time=`date +%s`
echo 4/5 paths redundancy cleaned in `expr $endclean_time - $startclean_time` seconds;



#COMPACTION OF SUPERREADS



startass_time=`date +%s`
../kMILL/src/kMILL short_read_connector_res.txt >>log 2>>log;

mv out_short_read_connector_res.txt.fa contigs2.fa
n50 contigs2.fa;


#~ rm noduplicate.fa notAligned.fa paths reads_corrected.fasta reads_corrected.unitigs.fa reads.fasta short_read_connector_res.txt >>log 2>>log;



