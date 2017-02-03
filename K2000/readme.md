K2000
=====

 

K2000 history
-------------

K2000 is derived from the KMill program written by [Camille
Marchet](http://people.rennes.inria.fr/Camille.Marchet/).

 

K2000 Purposes
--------------

K2000 is mainly used in the BWISE. It generates "phasitigs" from node ids of
read-coherent paths in the compacted de-Bruijn graph and the unitig sequences
corresponding to those node ids.

 

K2000 IOs
---------

### In

-   Set of read coherent paths in the compacted de-Bruijn graph. Each line is
    composed of id of compacted nodes (unitigs) seperated by ';'. A negative id
    designs the reverse complement of unitig. Example 14;-4;21 designs the path
    in a compacted de-Bruijn graph using nodes 14, the reverse complement of the
    node 4 and finally the node 21.

-   Unitigs corresponding to each id. This is a fasta sequence file. On sequence
    = one sequence of any length over the ACGT alphabet. First sequence is the
    unitig id '1', the second is '2', and so on. The headers of this fasta file
    are lost and ignored

-   An integer value k that is the kmer size as used in the de-Bruijn graph

-   Name of the output gfa file

-   [Name of the output fasta file]

### Out

-   A [GFA](https://github.com/GFA-spec/) graph formated. Each node corresponds
    to a simple path among the graph of unitigs.

-   If requiered, a fasta file in which each node of the gfa is output in this
    file in a fasta fashion

 

K2000 Requirement
-----------------

python \> 3.0

 

K2000 quick starting and test.
------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
./run_K2000.sh tests/dbg_paths.txt tests/unitigs.fa 240 my_first_gfa.gfa
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may run the `./test.sh` script. It runs the previous command, check results,
and remove created files.

 

Running K2000
-------------

use the *run_K2000.sh* script.

 

Directory Scripts 
------------------

The directory script contains diverse subtools, mainly useful for results
validation.

-   script validate compactions.

    -   usage: `python validate_compactions.py my_compacted_SR_file
        my_original_SR_file `(with SR meaning "Super Reads", indicated as
        `dbg_paths.txt `in the previous example)

        -   This validates that each of the SR in the original_SR_file is
            present in a compacted SR.

    -   Note that this script does not verify that all compacted super reads are
        correctly compacted. However, a bad compaction is detected during the
        generation of compacted sequences: if the k-1 overlap between two SR is
        not equal a WARNING message is printed (but the program does not stop).

        -   Example:

            $$WARNING unitigs ids 18600 and 9258 do not overlap, while they do
            in [18600, 9258, 200662 ...]$$

-   Counting unitigs coverages from SR file. It may be useful to get information
    about unitigs coverages. This is the purpose of the script
    split_count_ids.py.

    -   usage: `python split_count_ids.py my_original_SR_file > flat_unitigs`

        -   This simply spreads out the unitigs ids in the `flat_unitigs` file.
            Afterwards, this file may be counted / sorted as presented in the
            `howto_get_histograms.txt` file, including the usage of the
            `histo.R` R script.

             

             

 

### Authors

[pierre.peterlongo\@inria.fr](mailto:pierre.peterlongo@inria.fr )

[camille.marchet\@inria.fr](camille.marchet@inria.fr )

[antoine.limasset\@irisa.fr](antoine.limasset@irisa.fr )

 

 
