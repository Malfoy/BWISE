#! /usr/bin/env bash

EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


python3 ${EDIR}/check_python3_or_greater.py
if [ $? -ne 0 ] ; then
echo "*** K2000 needs python3 or greater to be ran *** "
echo "Please use python3 and re-run K2000"
exit 1
fi

function help {
    echo "THIS IS K2000 HELP:"
    echo "run_K2000.sh, deal with a super read file (from bgreat) with abundances to transform it as an assembly (.fa) and optionally a graph file (.gfa)"
    echo "Usage: ./run_K20000.sh [OPTIONS]"
    echo -e "\tMANDATORY:"
    echo -e "\t\t -i input file containing dbg paths as a list of unitig ids, eg. on line looks like \"-1;24;198;\""
    echo -e "\t\t -f: output fasta file"
    echo -e "\t\t -g: output gfa file"
    echo -e "\t\t -u: input fasta file containing unitig sequences. Note that unitig id 1 (as indicated in the paths input file) corresponds to the first unitig. "
    echo -e "\t\t -k: kmer size. "


    echo -e "\tK2000 OPTIONS:"
    echo -e "\t\t -c: Minimal conflict overlap. \n\t\t\t With c=0: K2000 is exact. \n\t\t\t With c>0: K2000 becomes greedy, in this case if a path A could be extended either by B or C and B and C are not collinear, then if the size of the overlap (A,B)>c and the size of the overlap (A,C)<c, then compact A-B but not A-C. If both overlaps are bigger than c or both smaller than c, no compaction is made. \n\t\t     Note that with c>0, size of unitigs has to be computable, thus K2000 needs to know the k value and the unitig length. Thus, with c>0, options -k and --unitig_file  are mandatory. [DEFAULT 0]"

    echo -e "\t\t -t: Dead end smaller or equal than this value are removed from the path graph.\n\t\t     Note that with C>0, size of unitigs has to be computable, thus K2000 needs to know the k value and the unitig length. Thus, with C>0, options -k and --unitig_file  are mandatory. [DEFAULT 0]"

    echo -e "\t\t -b: crush bubbles -- experimental"
    echo -e "\t\t -h: Prints this message and exist"
    echo "Any further question: read the readme file or contact us: bwise@inria.fr"
}



function run_full_k2000 { # Run the pipeline


    echo "*** REMOVE DUPLICATES AND COMPACT MAXIMAL SUPER READS *******"
    cmd="python3 ${EDIR}/K2000.py  -u ${unitig_file} -k ${k_value} -c ${min_conflict_overlap} -t ${max_tip}  ${dbg_path_file} -b ${crush}"
    echo ${cmd} " > ${original_dbg_path_file}_compacted_${min_conflict_overlap}"
    ${cmd} > ${original_dbg_path_file}_compacted_${min_conflict_overlap}
    if [ $? -ne 0 ] ; then
       (>&2 echo "There was a problem in the unitig compaction, K2000 ABORDED")
       exit 1
    fi



    echo "*** GENERATE GFA GRAPH FROM COMPACTED MAXIMAL SUPER READS ***"
    cmd="python3 ${EDIR}/K2000_msr_to_gfa.py ${original_dbg_path_file}_compacted_${min_conflict_overlap} ${unitig_file} ${k_value} ${dbg_path_file}"
    echo ${cmd} " > ${out_gfa}"
    ${cmd} > ${out_gfa}
    if [ $? -ne 0 ] ; then
           (>&2 echo "There was a problem in the unitig compaction during the GFA construction, K2000 ABORDED")
           exit 1
    fi


    echo "*** GENERATE FASTA FILE ***"
    cmd="python3 ${EDIR}/K2000_gfa_to_fasta.py ${out_gfa}"
    echo ${cmd} "> ${out_fasta}"
    $cmd > ${out_fasta}
    if [ $? -ne 0 ] ; then
           (>&2 echo "There was a problem in the unitig compaction during the Fasta construction, K2000 ABORDED")
    exit 1
    fi

}

# Initialize our own variables:

crush=0
min_conflict_overlap=0
max_tip=0
while getopts "hi:f:g:c:t:u:k:b" opt; do
    case "$opt" in
    h)
        help
        exit 0
    ;;
    f) out_fasta=$OPTARG
    ;;
    g) out_gfa=$OPTARG
    ;;
    i) dbg_path_file=$OPTARG
    ;;
    c) min_conflict_overlap=$OPTARG
    ;;
    t) max_tip=$OPTARG
    ;;
    u) unitig_file=$OPTARG
    ;;
    k) k_value=$OPTARG
    ;;
    b)
        crush=1
        ;;
    esac
done

if [ -z "$dbg_path_file" ]; then
    (>&2 echo "ERROR: You must provide a dbg_path_file (-i)")
    help
    exit
fi


if [ -z "$unitig_file" ]; then
    (>&2 echo "ERROR: You must provide a unitig_file (-u)")
    help
    exit
fi

if [ -z "$k_value" ]; then
    (>&2 echo "ERROR: You must provide a k value (-k)")
    help
    exit
fi

if [ -z "$out_fasta" ]; then
    (>&2 echo "ERROR: You must provide a out_fasta (-f)")
    help
    exit
fi

if [ -z "$out_gfa" ]; then
    (>&2 echo "ERROR: You must provide a out_gfa (-g)")
    help
    exit
fi





original_dbg_path_file=dbg_path_file
run_full_k2000

#~ if [ $crush -eq 1 ] ; then
    #~ ############### GENERATION OF CRUSHED BUBBLE GRAPH. UNSAFE AND COMMENTED FOR NOW ###############
    #~ echo "*** GENERATE A SIMPLIFIED CONSENSUS GRAPH (WARNING UNVALIDATED CODE)***"
    #~ cmd="python3 ${EDIR}/K2000_msr_to_simplified_msr.py ${original_dbg_path_file}_compacted_${min_conflict_overlap} ${unitig_file} ${k_value}"
    #~ echo ${cmd} "> ${original_dbg_path_file}_compacted_${min_conflict_overlap}_crushed"
    #~ ${cmd} > ${original_dbg_path_file}_compacted_${min_conflict_overlap}_crushed
    #~ if [ $? -ne 0 ] ; then
        #~ (>&2 echo "There was a problem in the graph simplification, K2000 ABORDED")
        #~ exit 1
    #~ fi
    #~ dbg_path_file=${original_dbg_path_file}_compacted_${min_conflict_overlap}_crushed
    #~ out_fasta=crushed_${out_fasta}
    #~ out_gfa=crushed_${out_gfa}
    #~ run_full_k2000

#~ fi

echo "*** K2000 DONE results are in files ${out_fasta} and ${out_gfa} ***"
exit 0
