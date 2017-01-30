#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
'''
Creation of a GFA file from a set of compacted maximal super reads 
@author  pierre peterlongo pierre.peterlongo@inria.fr
'''            

import sys
import getopt
import K2000_common as kc

def get_size_super_read_in_u(u,unitigs,k):
    '''compute the size of unitigs of a path stored in u (sum(size unitigs) - (|u|-1)*(k-1)'''
    cumulated_size_unitigs=0
    for unitig_id in u:
        if unitig_id<0: unitig_id=-unitig_id
        cumulated_size_unitigs+=len(unitigs[unitig_id-1])         # Ids starts at 1 (avoiding the -0=0 problem). Thus unitig i corresponds to unitigs[i-1]
    cumulated_size_unitigs-=(len(u)-1)*(k-1)                      # remove the size of the overlapping contigs
    return cumulated_size_unitigs

def show_right_edges (SR,x,id_x,unitigs,k,indexed_nodes):
    ''' Main function. For a given super read x, we find y that overlap x, and we print the links in a GFA style:
    L	11	+	12	-	overlap size
    Note that one treat x only if its canonical. 
    Four cases :
    1/ x overlaps y, with y canonical. One prints x + y + blabla
    2/ x overlaps y_, with y_ non canonical. In this case, y_ overlaps x. One of the two solutions has to be chosen. We chose min(idx,idy) (with idx,idy being the ids of the MSR x,y in SR) One searches the id of y, and one prints x + y - blabla. 
    3/ x_ overlaps y. same as 2. 
    4/ x_ overlaps y_. We do nothing, this case is treated when the entry of the function is y that thus overlaps x. 
    '''
    print (" x is ",x)
    if not kc.is_canonical(x): return
    n=len(x)
    
    
    # CASES 1 AND 2
    strandx='+'
    for len_u in range(1,n): # for each possible x suffix
        u=x[-len_u:]
        Y=SR.get_lists_starting_with_given_prefix(u)
        # if x in Y: Y.remove(x)   # we remove x itself from the list of y : note that it should not occur.
        # if x in Y: print ("ho ho ho",x)
        # print (x)
        # assert(x not in Y)
        if len(Y)==0: continue  # No y starting with u
        for y in Y:
            # detect the y strand 
            if kc.is_canonical(y): # CASE 1/
                strandy ='+'
                id_y=indexed_nodes.index(y)                                # get the id of the target node
            else:                  # CASE 2/
                strandy='-'
                id_y = indexed_nodes.index(kc.get_reverse_sr(y))
                if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y. 
            size_super_read = get_size_super_read_in_u(u,unitigs,k)
            # print the edges
            print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(size_super_read))
    
    
    # CASES 3 AND 4
    strandx='-'
    x_=kc.get_reverse_sr(x)
    for len_u in range(1,n): # for each possible x suffix
        u=x_[-len_u:]
        Y=SR.get_lists_starting_with_given_prefix(u)
        # assert(x_ not in Y)
        if len(Y)==0: continue  # No y starting with u
        for y in Y:
            if kc.is_canonical(y): # CASE 3
                strandy ='+'
                id_y=indexed_nodes.index(y)                                # get the id of the target node
                # we determine min(id_x,id_y)
                if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y. 
                size_super_read = get_size_super_read_in_u(u,unitigs,k)
                print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(size_super_read)) # note that strand x is always '-' and strandy is always '+' in this case. 
                
#            else: continue # CASE 4, nothing to do.

            
    
         
def print_GFA_edges(SR,unitigs,k,indexed_nodes):
    '''print each potiential edge in GFA format. Note that each edge is printed in one unique direction, the other is implicit'''
    for x_id in range(len(SR)):
        if x_id%100==0: sys.stderr.write("\t%.2f"%(100*x_id/len(SR))+"%\r")
        show_right_edges(SR,indexed_nodes[x_id],x_id,unitigs,k,indexed_nodes)
    sys.stderr.write("\t100.00%\n")
        
        
             
def print_GFA_nodes(SR, unitigs, k):
    '''print canonical unitigs'''
    node_id=-1
    for sr in SR.traverse():
        node_id+=1
        if node_id%100==0: sys.stderr.write("\t%.2f"%(100*node_id/len(SR))+"%\r")
        if kc.is_canonical(sr) :
            print ("S\t"+str(node_id)+"\t", end="")
            print_first_kmer=True
            for unitig_id in sr:
                reverse=False
                if unitig_id<0:                                         #the sequence is indicated as reverse complemented. Note that unitig ids start with 1, thus the -0 problem does not appear.
                    reverse=True
                    unitig_id=-unitig_id
                unitig=unitigs[unitig_id-1]                             # grab the good unitig. Ids starts at 1 (avoiding the -0=0 problem). Thus unitig i corresponds to unitigs[i-1]
                if reverse: unitig=kc.reverse_complement(unitig)        #reverse the untig if necessary
                if not print_first_kmer: unitig=unitig[k-1:]            #remove the k-1overlap
                if print_first_kmer: print_first_kmer=False             #do not print first kmer for next unitigs
                print (unitig,end="")
            print ()
    sys.stderr.write("\t100.00%\n")
            
            
def print_GFA_nodes_as_ids(SR, unitigs, k):
    '''print canonical unitigs ids'''
    node_id=-1
    for sr in SR.traverse():
        node_id+=1
        if is_canonical(sr) :
            print ("S\t"+str(node_id)+"\t", end="")
            for unitig_id in sr:
                print (str(unitig_id)+";", end="")
            print ()

def index_node_ids(SR):
    indexed_nodes=[]
    for sr in SR.traverse():
        indexed_nodes+=[sr]
    return indexed_nodes
        
    
             
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
    sys.stderr.write("Print GFA Nodes\n")
    print_GFA_nodes(SR,unitigs,k)
    sys.stderr.write("Print GFA Edges\n")
    indexed_nodes=index_node_ids(SR)
    print_GFA_edges(SR,unitigs,k, indexed_nodes)


if __name__ == "__main__":
     main()  
