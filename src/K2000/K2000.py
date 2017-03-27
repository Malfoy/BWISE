#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
'''
Compaction of super reads. Non ambiguous paths from a set of super reads are compacted
Resulting super reads are called MSR for Maximal Super Reads
@author  pierre peterlongo pierre.peterlongo@inria.fr
'''

import sys
import getopt
from multiprocessing.pool import ThreadPool as Pool
import K2000_common as kc
import sorted_list





def remove_y_strictly_included_in_x(x_ref,SR):
    ''' remove all y strictly included in x'''
    if len(x_ref)==1: return # as we removed strict equalities, no read can be included in a read of size one.
    n=len(x_ref)
    # print ("x",x)

    for x in [x_ref, kc.get_reverse_sr(x_ref)]:
        for position_suffix in range(0,n):
            u=x[position_suffix]
            Y=SR.get_lists_starting_with_given_prefix([u])
            # print (position_suffix, "-", x,"-",  u, "-" , Y)
            if x in Y: Y.remove(x)
            for y in Y:
                if len(y)+position_suffix <= n and kc.colinear(x,[y],[position_suffix]):
                    # print ("remove ",y, "in ",x,"pos",position_suffix)
                    SR.remove(y)
                    SR.remove(kc.get_reverse_sr(y))







def remove_strict_inclusions(SR):
    ''' remove all super reads strictly included in any other '''
    n = len(SR)
    # to_remove=[False for i in range(n)]
    checked=0

    for sr in SR.traverse():
        if checked%100==0: sys.stderr.write("      Removing inclusions, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\r")
        checked+=1
        # pool.apply_async(remove_y_strictly_included_in_x(sr,SR,to_remove))

        remove_y_strictly_included_in_x(sr,SR)#,to_remove)

    sys.stderr.write("      Removing inclusions, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\n")
    return SR




def get_len_ACGT(sr,unitig_lengths,k):
    ''' provides the cumulated length of the unitigs of a sr '''
    lenACGT = k
    for unitig_ids in sr:
        if unitig_ids>0:
            lenACGT+=unitig_lengths[unitig_ids-1]-k       # add the ACGT len of the corresponding unitig. -1 is due to the fact that unitigs start at zero while sr indices start at one.
        else:
            lenACGT+=unitig_lengths[-unitig_ids-1]-k

    return lenACGT



def right_unique_extention(SR,sr, unitig_lengths,k,gready):
    ''' return the unique  possible right sr extension with the largest overlap
        return None if no right extensions or if more than one possible non colinear extensions
    '''

    n=len(sr)
    #  **** Get the largest right overlap with sr ****
    for len_u in range(n-1,0,-1):
        u=sr[-len_u:]
        Y=SR.get_lists_starting_with_given_prefix(u)
        if sr in Y: Y.remove(sr)            # possible if x is repeated 2,2,2,2 for instance or if prefix = suffix (1,2,1 for instance)
        if len(Y)==0: continue              # No y starting with u
        if len(Y)>1: return None,len_u      # More than one unique y starting with u, for instance y and y'. Knowing that y is not included in y' it means necessary that y and y' are not colinear and that x is thus right extensible.
        y=Y[0]                              # We found the largest y right overlapping sr.

        lenACGT_u = get_len_ACGT(u,unitig_lengths,k)
        # **** check that all other right extensions are collinear with y.
        # get all y' such that LCSP(sr,y') in [1,len(u)-1]
        # get all y' starting with a suffix of u
        Y=[]
        starting_positions=[]
        for starting_suffix_position in range(1,len_u):
            suffix_u=u[starting_suffix_position:]
            # GREADY PART: don't check colinearity with sr that overlap "not enough x" with respect to largest possible overlap with x.
            # eg:
            # x:        ------------------
            # y:          --------------------
            # z:                        ------------------
            #             <-   >500   ->
            # in this case we don't check the y and z colinearity
            lenACGT_suffix_u = get_len_ACGT(suffix_u,unitig_lengths,k)          # TODO: optimize this.
            # print("diff is",lenACGT_u - lenACGT_suffix_u)
            if gready and lenACGT_u - lenACGT_suffix_u > 500:                              # TODO: this value should be a parameter and maybe a ratio.
                break
            # END OF THE GREADY PART
            others = SR.get_lists_starting_with_given_prefix(suffix_u)
            if len(others) >0:
                Y+=others
                starting_positions+=[starting_suffix_position for zz in range(len(others))]
        if len(starting_positions)>0 and not kc.colinear(y,Y,starting_positions): return None,len_u
        return y,len_u
    return None,None

def fusion (SR,x, unitig_lengths,k,gready):
    '''Main function. For a given super read x, we find y that overlap x with the highest overlap, such that :
    1/ there exists no other y' right overlapping x that is not collinear with y
    2/ there exists no other x' left overlapping y that is not collinear with x
    Once done, we compact x and y, and this is finished for x.
    '''

    y,len_u=right_unique_extention(SR,x, unitig_lengths,k,gready)               # Define, if exists, the unique y != x having the largest right overlap with x.
    if y==None: return 0                                                        # if no unique right extension, finished, x is not right extensible.
    y_= kc.get_reverse_sr(y)                                                    # what are left extentions of y, we right extend its reverse complement.
    xprime_, dontcare = right_unique_extention(SR,y_, unitig_lengths,k,gready)  # Define, if exists, the unique xprime_ (!= y_) having the largest right overlap with y_.
    if xprime_==None: return 0                                                  # if no unique left extension of the unique right extention of x, finished, x is not right extensible.


    if gready:                                                                  # if gready its possible* that y_ is extended with something else than x_
        if x!=kc.get_reverse_sr(xprime_): return 0                              # In this case, do not make the x/y fusion.
    else: 
        assert(x==kc.get_reverse_sr(xprime_))                                   # if not gready the unique right extension of reverse complement of x is x_

    # * eg:
    # x  -----------
    # y           ----------------
    # x'   /------------------         # x is extended only with y and y_ is extended only with x'_ (as the x/y overlap is too small
    # Note that latter: y will be compacted with x' (when the entry of the fusing function is y'_


    # ***** FUSION *****
    SR.remove(x)
    if not kc.is_palindromic(x):   SR.remove(xprime_)


    SR.remove(y)
    if not kc.is_palindromic(y):   SR.remove(kc.get_reverse_sr(y))

    new=x+y[len_u:]
    SR.sorted_add(new)
    if not kc.is_palindromic(new): SR.sorted_add(kc.get_reverse_sr(new))
    return 1


def compaction(SR, unitig_lengths,k,gready):
    ''' Compaction of all sr in SR
    If a sr was not compacted, i will never be compacted.
    If it was compacted, maybe it will be re-compacted later on. However, no need to re-run the fusion on this read as
     - either I could have been compacted on the left and this was done before or this will be done latter or
     - it will be right extended later: the 'new' (the compacted sr) sr is higher in the lexicographic order than the original sr (as it is longer), thus its position is higher in the SR data structure, thus it will be seen again later.
    Note that this may be not true in parallel computations.
    '''

    checked=0
    compacted=0
    n = len(SR)
    for sr in SR.traverse():
        if checked%100==0:
            sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\r")
        checked+=1
        witness = fusion(SR,sr,unitig_lengths,k,gready)
        if witness == 1: # a fusion was done
            compacted+=1
    sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\n")
    return SR


def main():
    '''
    Compaction of set of super reads coded as set of ids of unitigs
    '''

    k = int(sys.argv[3])

    gready = True
    if len(sys.argv)==5 and sys.argv[4]=="-e":
        gready = False
        sys.stderr.write("  Exact K2000 \n")
    else:
        sys.stderr.write("  Gready K2000 \n")

    sys.stderr.write("  Load unitig lengths \n")
    unitig_lengths = kc.load_unitig_lengths (sys.argv[2])

    sys.stderr.write("  Load super reads \r")
    SR=kc.generate_SR(sys.argv[1])
    sys.stderr.write("  Load super reads. Done - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Add reverse complements \r")
    kc.add_reverse_SR(SR)
    sys.stderr.write("  Add reverse complements. Done - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Conserve only unique super reads \r")
    SR.unique()
    sys.stderr.write("  Conserve only unique super reads. Done - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Sorting \r")
    SR.sort()
    sys.stderr.write("  Sorting. Done - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Remove strict inclusions\r")
    SR=remove_strict_inclusions(SR)
    sys.stderr.write("  Remove strict inclusions. Done - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Compaction of simple paths \r")
    SR=compaction(SR, unitig_lengths,k,gready)
    sys.stderr.write("  Compaction of simple paths. Done - nb SR="+ str(len(SR))+"\n")


    sys.stderr.write("  Print maximal super reads - nb canonical MSR="+ str(int(len(SR)/2))+"\n")
    kc.print_maximal_super_reads(SR)

if __name__ == "__main__":
     main()
