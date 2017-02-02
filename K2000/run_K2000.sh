if (( $# < 4 )); then
       echo "Need 4 or 5 parameters in the following order dbg_path_file unitig_file k_value out_file_gfa [out_file_fasta]"
       echo "see the readme file"
       exit
fi



EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )



in_sr=$1
in_unitigs=$2
in_k=$3
out_gfa=$4
echo "*** REMOVE DUPLICATES AND COMPACT MAXIMAL SUPER READS *******"
python ${EDIR}/K2000.py ${in_sr} > ${in_sr}_compacted
echo "*** GENERATE GFA GRAPH FROM COMPACTED MAXIMAL SUPER READS ***"
python ${EDIR}/K2000_msr_to_gfa.py ${in_sr}_compacted ${in_unitigs} ${in_k} #> ${out_gfa}
if (( $# == 5 )); then

       echo "*** GENERATE FASTA FILE ***"
       out_fasta=$5
       python ${EDIR}/K2000_gfa_to_fasta.py ${out_gfa} > ${out_fasta}
fi

echo "*** K2000 DONE ***"