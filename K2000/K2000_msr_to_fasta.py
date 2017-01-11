#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

'''
From a GFA file representing a graph of phasitigs, output only node sequences in a fasta fashion.
@author  pierre peterlongo pierre.peterlongo@inria.fr
'''                

import sys
import getopt
import K2000_common as kc
   
        
             
def print_GFA_nodes_as_fasta(SR, unitigs, k):
    '''print canonical unitigs'''
    node_id=-1
    for sr in SR:
        node_id+=1
        if kc.is_canonical(sr) :
            print (">node_"+str(node_id))
            print_first_kmer=True
            for unitig_id in sr:
                reverse=False
                if unitig_id<0:                                         #the sequence is indicated as reverse complemented. Note that unitig ids start with 1, thus the -0 problem does not appear.
                    reverse=True
                    unitig_id=-unitig_id
                unitig=unitigs[unitig_id-1]                             # grab the good unitig. Ids starts at 1 (avoiding the -0=0 problem). Thus unitig i corresponds to unitigs[i-1]
                if reverse: unitig=kc.reverse_complement(unitig)           #reverse the untig if necessary
                if not print_first_kmer: unitig=unitig[k-1:]            #remove the k-1overlap
                if print_first_kmer: print_first_kmer=False             #do not print first kmer for next unitigs
                print (unitig,end="")
            print ()
            
    
            

def main():
    '''
    Creation of a GFA file from a set of compacted maximal super reads 
    '''
    # SR=[[1,3,4],[14],[4,6],[-6,-1,-2],[4,6,7,9],[9,10,11,12],[4,6],[-13, -12, -11]]
    
    SR=kc.generate_SR(sys.argv[1])
    unitigs=kc.load_unitigs(sys.argv[2])
    k = int(sys.argv[3])
    kc.add_reverse_SR(SR)
    SR.sort()
    
    print_GFA_nodes_as_fasta(SR,unitigs,k)


if __name__ == "__main__":
     main()  
