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

def get_sr_id(sr):
    ''' returns the id of a sr 
    WARNING: here sr contains as last value its unique id. 
    '''
    return int(sr[-1].split("_")[1])
    
def get_reverse_sr_id(sr,SR):
    ''' returns the id of the reverse complement of sr 
    1/ find the sequence of the sr_ in SR
    2/ grab its id
    WARNING: here sr contains as last value its unique id. 
    '''
    #1/
    without_id_reverse_sr = kc.get_reverse_sr(sr[:-1])                  # get the sr reverse complement
    Y=SR.get_lists_starting_with_given_prefix(without_id_reverse_sr)    # find the reverse complement in the list. 
    for y in Y:                                                         # should be of size >= 1. One of them is exactly 'without_id_reverse_sr' plus its id. 
        if len(y) == len(without_id_reverse_sr)+1:                      # 'y' is 'without_id_reverse_sr' with its node id
            return get_sr_id(y)                                         # 2/
    return None                                                         # Should not happend

def get_size_super_read_in_u(u,unitigs,k):
    '''compute the size of unitigs of a path stored in u (sum(size unitigs) - (|u|-1)*(k-1)'''
    cumulated_size_unitigs=0
    for unitig_id in u:
        if unitig_id<0: unitig_id=-unitig_id
        cumulated_size_unitigs+=len(unitigs[unitig_id-1])         # Ids starts at 1 (avoiding the -0=0 problem). Thus unitig i corresponds to unitigs[i-1]
    cumulated_size_unitigs-=(len(u)-1)*(k-1)                      # remove the size of the overlapping contigs
    return cumulated_size_unitigs

def show_right_edges (SR,x,id_x,unitigs,k):
    ''' Main function. For a given super read x, we find y that overlap x, and we print the links in a GFA style:
    L	11	+	12	-	overlap size
    Note that one treat x only if its canonical. 
    Four cases :
    1/ x overlaps y, with y canonical. One prints x + y + blabla
    2/ x overlaps y_, with y_ non canonical. In this case, y_ overlaps x. One of the two solutions has to be chosen. We chose min(idx,idy) (with idx,idy being the ids of the MSR x,y in SR) One searches the id of y, and one prints x + y - blabla. 
    3/ x_ overlaps y. same as 2. 
    4/ x_ overlaps y_. We do nothing, this case is treated when the entry of the function is y that thus overlaps x. 
    WARNING: here x and each sr in SR contain as last value its unique id. 
    '''
    x=x[:-1]                                # remove the x_id from the x sr
    if not kc.is_canonical(x): return
    n=len(x)
    
    # CASES 1 AND 2
    strandx='+'
    # print ("x is", x)
    for len_u in range(1,n): # for each possible x suffix
        u=x[-len_u:]
        # print ("u is", u)
        Y=SR.get_lists_starting_with_given_prefix(u)
        # if x in Y: Y.remove(x)            # we remove x itself from the list of y : note that it should not occur.
        # print (x)
        # assert(x not in Y)
        if len(Y)==0: continue              # No y starting with u
        for y in Y:
            # detect the y strand 
            # CASE 1/
            if kc.is_canonical(y[:-1]):     # remove the last value that corresponds to the node id
                strandy ='+'
                # id_y=indexed_nodes.index(y)                               # get the id of the target node
                id_y=get_sr_id(y)                                           # last value is the node id, here the id of the target node
            # CASE 2/
            else:                       
                strandy='-'
#                id_y = indexed_nodes.index(kc.get_reverse_sr(y))
                id_y=get_reverse_sr_id(y,SR)                                # find the reverse of list y in SR to grab its id. 
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
            if kc.is_canonical(y[:-1]): # CASE 3
                strandy ='+'
               # id_y=indexed_nodes.index(y)                                # get the id of the target node
                id_y=get_sr_id(y)                                           # last value is the node id, here the id of the target node
                # we determine min(id_x,id_y)
                if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y. 
                size_super_read = get_size_super_read_in_u(u,unitigs,k)
                print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(size_super_read)+"M") # note that strand x is always '-' and strandy is always '+' in this case. 
                
#            else: continue # CASE 4, nothing to do.

            
    
         
def print_GFA_edges(SR,unitigs,k):
    '''print each potiential edge in GFA format. Note that each edge is printed in one unique direction, the other is implicit
    WARNING: here each sr in SR contains as last value its unique id. 
    '''
    for sr in SR.traverse():
        x_id = get_sr_id(sr)                                         # last value is the node id
        if x_id%100==0: sys.stderr.write("\t%.2f"%(100*x_id/len(SR))+"%\r")
        show_right_edges(SR,sr,x_id,unitigs,k)
    sys.stderr.write("\t100.00%\n")
        
        
             
def print_GFA_nodes(SR, unitigs, size_overlap,unitig_coverages):
    '''print canonical unitigs
    WARNING: here each sr in SR contains as last value its unique id. 
    '''
    nb_errors=0
    for sr in SR.traverse():
        node_id = get_sr_id(sr)                                         # last value is the node id
        sr = sr[:-1]                                                    # remove the last value that corresponds to the node id
        if node_id%100==0: sys.stderr.write("\t%.2f"%(100*node_id/len(SR))+"%\r")
        if kc.is_canonical(sr):              
            print ("S\t"+str(node_id)+"\t", end="")
            print_first_kmer=True
            previous_id=-1
            previous_overlap=""                                         #used only to assert a good k-1 overlap. 
            sum_coverage=0
            for unitig_id in sr:
                reverse=False
                # print ("\n",str(unitig_id-1))
                if unitig_id<0:                                         #the sequence is indicated as reverse complemented. Note that unitig ids start with 1, thus the -0 problem does not appear.
                    reverse=True
                    unitig_id=-unitig_id
                unitig=unitigs[unitig_id-1]                             # grab the good unitig. Ids starts at 1 (avoiding the -0=0 problem). Thus unitig i corresponds to unitigs[i-1]
                sum_coverage+=(unitig_coverages[unitig_id]*len(unitig)) # The unitig had been seen unitig_coverages[unitig_id] times. As downstreat visualization tool as bandage divides by the sequence size, we multiply the coverage to the unitig by the size of the unitig
                if reverse: unitig=kc.reverse_complement(unitig)        #reverse the untig if necessary
                if previous_overlap != "":                              # overlap validation
                    if(unitig[:size_overlap] != previous_overlap):      # overlap validation
                        nb_errors+=1
                        # sys.stderr.write("\n WARNING unitigs ids "+ str(previous_id)+" and "+str(unitig_id)+" do not overlap, while they do in "+str(sr)+"\n")
                        # sys.stderr.write("Suffix "+previous_overlap+" size="+str(len(previous_overlap))+"\n")
                        # sys.stderr.write("Prefix "+unitig[:size_overlap]+"\n")
                previous_overlap = unitig[-size_overlap:]               # store the suffix of size k-1 to check the next overla
                previous_id = unitig_id
                if not print_first_kmer: unitig=unitig[size_overlap:]   # remove the k-1 overlap
                print_first_kmer=False                                  #do not print first kmer for next unitigs
                print (unitig,end="")
            print ("\tFC:i:"+str(sum_coverage),end="")                         # Fragment count
            # print comments about untigs 
            print ("\t# Unitigs/coverage: ",end="")
            for unitig_id in sr:
                if unitig_id<0: unitig_id=-unitig_id
                print (str(unitig_id)+"/"+str(unitig_coverages[unitig_id])+" ",end="")
            print ()
                
    sys.stderr.write("\t100.00% -- "+ str(nb_errors)+" error(s)\n" )
            
            
def print_GFA_nodes_as_ids(SR, unitigs, k, unitig_coverages):
    '''print canonical unitigs ids
    WARNING: here each sr in SR contains as last value its unique id. 
    '''
    for sr in SR.traverse():
        node_id = get_sr_id(sr)                        # last value is the node id
        sr = sr[:-1]                                   # remove the last value that corresponds to the node id
        if  kc.is_canonical(sr):                       # remove the last value that corresponds to the node id
            print ("S\t"+str(node_id)+"\t", end="")
            for unitig_id in sr:                       # remove the last value that corresponds to the node id
                print (str(unitig_id)+";", end="")
            print ()


def load_unitig_coverage(file_name):
    sr_file = open(file_name, 'r')
    unitig_coverages={}
    for line in sr_file:
        
        line = line.rstrip()[:-1].split(';')
        for unitig_id in line: 
            unitig_id=int(unitig_id)
            if unitig_id<0: unitig_id=-unitig_id
            if unitig_id not in unitig_coverages: unitig_coverages[unitig_id]=0
            unitig_coverages[unitig_id]+=1
    return unitig_coverages

def main():
    '''
    Creation of a GFA file from a set of compacted maximal super reads 
    '''
    # SR=[[1,3,4],[14],[4,6],[-6,-1,-2],[4,6,7,9],[9,10,11,12],[4,6],[-13, -12, -11]]
    
    SR=kc.generate_SR(sys.argv[1])
    unitigs=kc.load_unitigs(sys.argv[2])
    k = int(sys.argv[3])
    unitig_coverages=load_unitig_coverage(sys.argv[4])
    kc.add_reverse_SR(SR)
    SR.sort()
    SR.index_nodes()                          # This adds a final value to each sr, providing its node id. 
    sys.stderr.write("Print GFA Nodes\n")
    print_GFA_nodes(SR,unitigs,k-1,unitig_coverages)
    # print_GFA_nodes_as_ids(SR,unitigs,k-1)
    sys.stderr.write("Print GFA Edges\n")
    print_GFA_edges(SR,unitigs,k)


if __name__ == "__main__":
     main()  
