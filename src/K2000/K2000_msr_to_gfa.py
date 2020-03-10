#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
Creation of a GFA file from a set of compacted maximal super reads
@author  pierre peterlongo pierre.peterlongo@inria.fr
'''

import sys
import getopt
import K2000_common as kc

# def get_size_super_read_in_u(u,unitigs,k):
#     '''compute the size of unitigs of a path stored in u (sum(size unitigs) - (|u|-1)*(k-1)'''
#     cumulated_size_unitigs=0
#     for unitig_id in u:
#         if unitig_id<0: unitig_id=-unitig_id
#         cumulated_size_unitigs+=len(unitigs[unitig_id-1])         # Ids starts at 1 (avoiding the -0=0 problem). Thus unitig i corresponds to unitigs[i-1]
#     cumulated_size_unitigs-=(len(u)-1)*(k-1)                      # remove the size of the overlapping contigs
#     return cumulated_size_unitigs

# def show_right_edges (MSR,x,id_x,unitigs,k):
#     ''' Main function. For a given super read x, we find y that overlap x, and we print the links in a GFA style:
#     L    11    +    12    -    overlap size
#     Note that one treat x only if its canonical.
#     Four cases :
#     1/ x overlaps y, with y canonical. One prints x + y + blabla
#     2/ x overlaps y_, with y_ non canonical. In this case, y_ overlaps x. One of the two solutions has to be chosen. We chose min(idx,idy) (with idx,idy being the ids of the MSR x,y in SR) One searches the id of y, and one prints x + y - blabla.
#     3/ x_ overlaps y. same as 2.
#     4/ x_ overlaps y_. We do nothing, this case is treated when the entry of the function is y that thus overlaps x.
#     WARNING: here x and each msr in MSR contain as last value its unique id.
#     '''
#     x=x[:-1]                                # remove the x_id from the x msr
#     if not kc.is_canonical(x): return
#     n=len(x)
#
#     # CASES 1 AND 2
#     strandx='+'
#     # print ("x is", x)
#     for len_u in range(1,n): # for each possible x suffix
#         u=x[-len_u:]
#         # print ("u is", u)
#         Y=MSR.get_lists_starting_with_given_prefix(u)
#         # if x in Y: Y.remove(x)            # we remove x itself from the list of y : note that it should not occur.
#         # print (x)
#         # assert(x not in Y)
#         if len(Y)==0: continue              # No y starting with u
#         for y in Y:
#             # detect the y strand
#             # CASE 1/
#             if kc.is_canonical(y[:-1]):     # remove the last value that corresponds to the node id
#                 strandy ='+'
#                 # id_y=indexed_nodes.index(y)                               # get the id of the target node
#                 id_y=kc.get_msr_id(y)                                           # last value is the node id, here the id of the target node
#             # CASE 2/
#             else:
#                 strandy='-'
# #                id_y = indexed_nodes.index(kc.get_reverse_msr(y))
#                 id_y=kc.get_reverse_msr_id(y,MSR)                                # find the reverse of list y in MSR to grab its id.
#                 if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y.
#             size_super_read = get_size_super_read_in_u(u,unitigs,k)
#             # print the edges
#             print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(size_super_read)+"M")
#
#
#     # CASES 3 AND 4
#     strandx='-'
#     x_=kc.get_reverse_sr(x)
#     for len_u in range(1,n): # for each possible x suffix
#         u=x_[-len_u:]
#         Y=MSR.get_lists_starting_with_given_prefix(u)
#         # assert(x_ not in Y)
#         if len(Y)==0: continue  # No y starting with u
#         for y in Y:
#             if kc.is_canonical(y[:-1]): # CASE 3
#                 strandy ='+'
#                # id_y=indexed_nodes.index(y)                                # get the id of the target node
#                 id_y=kc.get_msr_id(y)                                           # last value is the node id, here the id of the target node
#                 # we determine min(id_x,id_y)
#                 if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y.
#                 size_super_read = get_size_super_read_in_u(u,unitigs,k)
#                 print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(size_super_read)+"M") # note that strand x is always '-' and strandy is always '+' in this case.
#
# #            else: continue # CASE 4, nothing to do.




# def print_GFA_edges(MSR,unitigs,k):
#     '''print each potiential edge in GFA format. Note that each edge is printed in one unique direction, the other is implicit
#     WARNING: here each msr in MSR contains as last value its unique id.
#     '''
#     for msr in MSR.traverse():
#         x_id = kc.get_msr_id(msr)                                         # last value is the node id
#         if x_id%100==0: sys.stderr.write("\t%.2f"%(100*x_id/len(MSR))+"%\r")
#         show_right_edges(MSR,msr,x_id,unitigs,k)
#     sys.stderr.write("\t100.00%\n")



# def print_GFA_nodes(MSR, unitigs, size_overlap,number_mapped_sr):
#     '''print canonical unitigs
#     WARNING: here each msr in MSR contains as last value its unique id.
#     '''
#     nb_errors=0
#     for indexed_msr in MSR.traverse():
#         node_id = kc.get_msr_id(indexed_msr)                              # last value is the node id
#         msr = indexed_msr[:-1]                                            # remove the last value that corresponds to the node id
#         if node_id%100==0: sys.stderr.write("\t%.2f"%(100*node_id/len(MSR))+"%\r")
#         if not kc.is_canonical(msr): continue
#
#         print ("S\t"+str(node_id)+"\t", end="")
#         print_first_kmer=True
#         previous_id=-1
#         previous_overlap=""                                         #used only to assert a good k-1 overlap.
#         size_msrACGT=0
#         msrACGT=""
#         for unitig_id in msr:
#             reverse=False
#             # print ("\n",str(unitig_id-1))
#             if unitig_id<0:                                         #the sequence is indicated as reverse complemented. Note that unitig ids start with 1, thus the -0 problem does not appear.
#                 reverse=True
#                 unitig_id=-unitig_id
#             unitig=unitigs[unitig_id-1]                             # grab the good unitig. Ids starts at 1 (avoiding the -0=0 problem). Thus unitig i corresponds to unitigs[i-1]
#             # sum_coverage+=(unitig_coverages[unitig_id]*len(unitig)) # The unitig had been seen unitig_coverages[unitig_id] times. As downstreat visualization tool as bandage divides by the sequence size, we multiply the coverage to the unitig by the size of the unitig
#             if reverse: unitig=kc.reverse_complement(unitig)        #reverse the untig if necessary
#             if previous_overlap != "":                              # overlap validation
#                 if(unitig[:size_overlap] != previous_overlap):      # overlap validation
#                     nb_errors+=1
#                     # sys.stderr.write("\n WARNING unitigs ids "+ str(previous_id)+" and "+str(unitig_id)+" do not overlap, while they do in "+str(msr)+"\n")
#                     # sys.stderr.write("Suffix "+previous_overlap+" size="+str(len(previous_overlap))+"\n")
#                     # sys.stderr.write("Prefix "+unitig[:size_overlap]+"\n")
#             previous_overlap = unitig[-size_overlap:]               # store the suffix of size k-1 to check the next overlap
#             previous_id = unitig_id
#             if not print_first_kmer: unitig=unitig[size_overlap:]   # remove the k-1 overlap
#             print_first_kmer=False                                  #do not print first kmer for next unitigs
#             # print (unitig,end="")
#             msrACGT+=unitig
#         # size_msrACGT=len(msrACGT)
#         # coverage=number_mapped_sr[node_id]/float(size_msrACGT)
#         print (msrACGT+"\tFC:i:"+str(number_mapped_sr[node_id]))#+"\t#AVG:f:%.2f"%(coverage))
#
#
#
#     sys.stderr.write("\t100.00% -- "+ str(nb_errors)+" error(s)\n" )


def print_GFA_nodes_as_ids(MSR, unitigs, k):
    '''print canonical unitigs ids
    WARNING: here each msr in MSR contains as last value its unique id.
    '''
    for msr in MSR.traverse():
        node_id = kc.get_msr_id(msr)                        # last value is the node id
        msr = msr[:-1]                                   # remove the last value that corresponds to the node id
        if  kc.is_canonical(msr):                       # remove the last value that corresponds to the node id
            print ("S\t"+str(node_id)+"\t", end="")
            for unitig_id in msr:                       # remove the last value that corresponds to the node id
                print (str(unitig_id)+";", end="")
            print ()

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def get_number_of_sr_occurring_in_a_msr(msr,SR,unitigs,size_overlap):
    '''
    given a msr,
    1/ find sr mapping on it (forward msr strand, but all sr are both forward and reverse so no need to compute reverse msr)
    2/ count the length of each sr
    3/ return the cumulative length of mapped sr.
    '''
    n=len(msr)
    nb_mapped=0

    for position_suffix in range(0,n):                                              #each possible possible unitig id from msr
        u=msr[position_suffix]
        anchored_sr_set=SR.get_lists_starting_with_given_prefix([u])                #get sr anchored to msr
        for sr in anchored_sr_set:                                                  # check if its colinear
            if len(sr)+position_suffix <= n and kc.colinear(msr,[sr],[position_suffix]):
                # len_this_sr=0                                                       # compute the len of the mapped sr
                # for unitig_id in sr:                                                # compute the len of the mapped sr
                    # unitig_id = int(unitig_id)                                      # compute the len of the mapped sr
                    # if unitig_id<0: unitig_id=-unitig_id                            # compute the len of the mapped sr
                    # len_this_sr += len(unitigs[unitig_id-1])                        # compute the len of the mapped sr. -1 because of indexing is zero based while numbers are one based.
                # len_this_sr -= size_overlap*(len(sr) -1)                            # remove ovelap size
                nb_mapped+=1                                                       # a new sr mapped
                # sys.stderr.write("\n"+str(msr)+" "+str(sr)+" "+str(len_this_sr)+" "+str(mapped_len))
    return nb_mapped



def compute_number_mapped(indexedMSR,SRfile_name,unitigs,size_overlap):
    ''' compute the cumulative mapped_len of each MSR
    1/ for each canonical msr
    2/ map all raw (potentially redundant and potentially strictly included) sr on msr (seen as forward or reverse)
    3/ count and return the cumulative length of mapped sr
    '''
    n=len(indexedMSR)
    SR=kc.generate_SR(SRfile_name)                                          #load and sort original sr file
    kc.add_reverse_SR(SR)                                                   #load and sort original sr file
    SR.sort()                                                               #load and sort original sr file
    number_mapped_sr={}                                                     #function results : cumulative len of mapped sr on each msr
    checked=0                                                               #cosmetics

    for msr in indexedMSR.traverse():                                       # 1/
        if checked%100==0: sys.stderr.write("\t%.2f"%(100*checked/n)+"%\r")
        checked+=1
        if not kc.is_canonical(msr[:-1]): continue
        msr_id = kc.get_msr_id(msr)
        number_mapped_sr[msr_id]= get_number_of_sr_occurring_in_a_msr(msr[:-1]                   ,SR,unitigs,size_overlap) # 2 and 3


    sys.stderr.write("\t%.2f"%(100*checked/n)+"%\n")
    return number_mapped_sr




def main():
    '''
    Creation of a GFA file from a set of compacted maximal super reads
    '''
    # SR=[[1,3,4],[14],[4,6],[-6,-1,-2],[4,6,7,9],[9,10,11,12],[4,6],[-13, -12, -11]]
    csr=kc.generate_SR(sys.argv[1])
    unitigs=kc.load_unitigs(sys.argv[2])
    k = int(sys.argv[3])



    csr.add_reverses()
    
    csr.sort()
    csr.index_nodes()                          # This adds a final value to each sr, providing its node id.
    # sys.stderr.write("Compute GFA Node coverages\n")
    # number_mapped_sr=compute_number_mapped(MSR,sys.argv[4],unitigs,k-1)
    # sys.stderr.write("Print GFA Nodes\n")
    csr.print_GFA(unitigs,k-1)#,number_mapped_sr)
    # print_GFA_nodes_as_ids(MSR,unitigs,k-1)
    # sys.stderr.write("Print GFA Edges\n")
    # print_GFA_edges(MSR,unitigs,k)


if __name__ == "__main__":
     main()
