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

### Out

A [GFA](https://github.com/GFA-spec/) graph formated. Each node corresponds to a
simple path among the graph of unitigs.

 

K2000 Requirement
-----------------

python \> 3.0

 

K2000 quick starting and test. 
-------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
./run_K2000.sh tests/dbg_paths.txt tests/unitigs.fa 220 my_first_gfa.gfa
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may run the `./test.sh` script. It runs the previous command, check results,
and remove created files.

 

Running K2000
-------------

use the *run_K2000.sh *script.

 

### Authors

pierre.peterlongo\@inria.fr

camille.marchet\@inria.fr

antoine.limasset\@irisa.fr

 

 
