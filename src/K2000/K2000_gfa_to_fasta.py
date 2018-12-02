#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

'''
From a GFA file representing a graph of phasitigs, output only node sequences in a fasta fashion.
@author  pierre peterlongo pierre.peterlongo@inria.fr
'''

import sys
import getopt
import K2000_common as kc



def print_GFA_nodes_as_fasta(gfa_file):
    '''print gfa nodes as fasta'''
    GFAs = open(gfa_file,"r")
    for line in GFAs:
        line = line.strip().split() #S       3       TTCGATAAATTGATCCAGGCTGCCGTCCAGCACGGCC...
        if line[0] != 'S': continue
        print(">node_"+line[1])
        print(line[2])

    GFAs.close()




def main():
    '''
    Creation of a fasta file from a GFA file.
    '''
    print_GFA_nodes_as_fasta(sys.argv[1])


if __name__ == "__main__":
     main()
