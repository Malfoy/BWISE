#! /usr/bin/env bash
if (( $# < 5 )); then
       echo "Need 5 parameters in the following order < dbg_path_file - unitig_file - k_value - out_file_gfa - out_file_fasta >"
       echo "see the readme file"
       exit
fi


EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


python3 ${EDIR}/check_python3_or_greater.py
if [ $? -ne 0 ] ; then
echo "*** K2000 needs python3 or greater to be ran *** "
echo "Please use python3 and re-run K2000"
exit 1
fi

in_sr=$1
original_in_sr=$1
in_unitigs=$2
in_k=$3
out_gfa=$4
echo "*** REMOVE DUPLICATES AND COMPACT MAXIMAL SUPER READS *******"
# BACKTOBWISE: UNCOMMENT THESE NEXT FIVE LINES
# python3 ${EDIR}/K2000.py ${in_sr} ${in_unitigs} ${in_k} -e > ${in_sr}_compacted
# if [ $? -ne 0 ] ; then
#        echo "There was a problem in the unitig compaction, K2000 ABORDED"
#        exit 1
# fi

# BACKTOBWISE: COMMENT THESE NEXT NINE LINES
for min_overlap in 51 101 151 201; do       # TODO: input list
    echo ${min_overlap}; 
    python3 ${EDIR}/K2000.py ${in_sr} ${in_unitigs} ${in_k} ${min_overlap} > ${original_in_sr}_compacted_${min_overlap}
    if [ $? -ne 0 ] ; then
       echo "There was a problem in the unitig compaction, K2000 ABORDED"
       exit 1
    fi
    in_sr=${original_in_sr}_compacted_${min_overlap}
done


echo "*** GENERATE GFA GRAPH FROM COMPACTED MAXIMAL SUPER READS ***"
python3 ${EDIR}/K2000_msr_to_gfa.py ${in_sr} ${in_unitigs} ${in_k} ${in_sr} > ${out_gfa}
if [ $? -ne 0 ] ; then
       echo "There was a problem in the unitig compaction during the GFA construction, K2000 ABORDED"
       exit 1
fi

echo "*** GENERATE FASTA FILE ***"
out_fasta=$5
python3 ${EDIR}/K2000_gfa_to_fasta.py ${out_gfa} > ${out_fasta}
if [ $? -ne 0 ] ; then
       echo "There was a problem in the unitig compaction during the Fasta construction, K2000 ABORDED"
exit 1
fi


echo "*** K2000 DONE ***"
exit 0
