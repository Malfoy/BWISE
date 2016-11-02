#!/bin/bash



readFileName=""
kmersize=62



while getopts "hr:" opt; do
case $opt in

h)
help
exit
;;


r)
echo "use reference file: $OPTARG" >&2
readFileName=$OPTARG
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

if [ -z "$readFileName" ]; then
help
exit 0
fi



#cleaning
mkdir folderAssembly;
cd folderAssembly;
rm *>>log 2>>log;



start_time=`date +%s`



#correct the sequences
startcorrection_time=`date +%s`
../bloocoo//build/bin/Bloocoo -file ../$readFileName  -kmer-size 31 -abundance-min 5 -out reads_corrected.fa >>log 2>>log;
#~ ../convertOneLineFasta.py reads_corrected.fasta > reads_corrected.fa
endcorrection_time=`date +%s`
echo 1/5 reads corrected in `expr $endcorrection_time - $startcorrection_time` seconds;



#dbg construction
startgraph_time=`date +%s`
../bcalm/build/bcalm -in reads_corrected.fa -kmer-size $kmersize -abundance-min 10 >>log 2>>log;
../kMILL/src/kMILL reads_corrected.unitigs.fa $((kmersize-1)) $((kmersize-2)) >>log 2>>log;
mv out_k$((kmersize-1))_.fa out.fa;
endgraph_time=`date +%s`
echo 2/5 graph constructed in `expr $endgraph_time - $startgraph_time` seconds;



#mapping of read on the dbg
startmap_time=`date +%s`
../BGREAT2/bgreat -k $kmersize -r reads_corrected.fa -g out.fa -t 20 -G -c -m 0 -e 1 >>log 2>>log;
endmap_time=`date +%s`
echo 3/5 read mapped on graph in `expr $endmap_time - $startmap_time` seconds;



#duplicate path elimination
startclean_time=`date +%s`
../kMILL/src/pathsCleaner paths 2	 >>log 2>>log;
#cat out.fa >> noduplicate.fa;



#elimination of contained path (non maximal one)
echo "noduplicate.fa" > bank
../BREADY/short_read_connector.sh -b bank -q bank >> log 2>>log;
rm *.h5;
endclean_time=`date +%s`
echo 4/5 paths redundancy cleaned in `expr $endclean_time - $startclean_time` seconds;



#compaction of maximal paths
startass_time=`date +%s`
../kMILL/src/kMILL short_read_connector_res.txt 2000 >>log 2>>log;
endass_time=`date +%s`
echo 5/5 maximal paths compacted in `expr $endass_time - $startass_time` seconds;
mv out_k2000_.fa contigs1.fa;



end_time=`date +%s`
echo Total execution time was `expr $end_time - $start_time` seconds;



#number of contigs obtained
#~ grep ">" -c cleanedContigs.fa;
#size of the contig file
../Count.py contigs1.fa;

rm noduplicate.fa notAligned.fa paths reads_corrected.fasta reads_corrected.unitigs.fa reads.fasta short_read_connector_res.txt >>log 2>>log;



