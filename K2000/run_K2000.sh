if (( $# != 4 )); then
    echo "Need 4 parameters in the following order dbg_path_file unitig_file k_value out_file "
    echo "see the readme file"
    exit
fi

in_sr=$1
in_unitigs=$2
in_k=$3
out_gfa=$4
python K2000.py ${in_sr} > compacted_unitigs.txt
python K2000_msr_to_gfa.py compacted_unitigs.txt ${in_unitigs} ${in_k} > ${out_gfa}