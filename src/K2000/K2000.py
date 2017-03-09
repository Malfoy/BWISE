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
                if len(y)+position_suffix <= n and colinear(x,[y],[position_suffix]):
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
    
    


    
def colinear(x,X,starting_positions):
    ''' Check that all sr in X are colinear with x
    For each sr in X, one knows its starting position on x, with table starting_positions'''
    for i in range(len(X)):
        other = X[i]
        starting_position = starting_positions[i]
        pos=0
        while True:
            if pos>=len(other) or pos+starting_position>=len(x) : break
            if other[pos] != x[pos+starting_position]:          # "non colinear"
                return False
            pos+=1
            
             
    return True

def right_unique_extention(SR,sr):
    ''' return the unique  possible right sr extension with the largest overlap
        return None if no right extensions or if more than one possible non colinear extensions
    '''

    n=len(sr)
    #  **** Get the largest right overlap with sr ****
    for len_u in range(n-1,0,-1): 
        u=sr[-len_u:]
        Y=SR.get_lists_starting_with_given_prefix(u)
        if sr in Y: Y.remove(sr)            # possible if x is repeated 2,2,2,2 for instance.
        if len(Y)==0: continue              # No y starting with u
        if len(Y)>1: return None,len_u      # More than one unique y starting with u, for instance y and y'. Knowing that y is not included in y' it means necessary that y and y' are not colinear and that x is thus right extensible. 
        y=Y[0]                              # We found the largest y right overlapping sr.
    
        # **** check that all other right extensions are collinear with y.
        # get all y' such that LCSP(sr,y') in [1,len(u)-1]
        # get all y' starting with a suffix of u
        Y=[]
        starting_positions=[]
        for starting_suffix_position in range(1,len_u):
            suffix_u=u[starting_suffix_position:]
            others = SR.get_lists_starting_with_given_prefix(suffix_u)
            if len(others) >0:
                Y+=others
                starting_positions+=[starting_suffix_position for zz in range(len(others))]
        if len(starting_positions)>0 and not colinear(y,Y,starting_positions): return None,len_u
        return y,len_u
    return None,None
    
def fusion (SR,x):
    '''Main function. For a given super read x, we find y that overlap x with the highest overlap, such that : 
    1/ there exists no other y' right overlapping x that is not collinear with y
    2/ there exists no other x' left overlapping y that is not collinear with x
    Once done, we compact x and y, and this is finished for x. 
    '''
    
    y,len_u=right_unique_extention(SR,x)                    # Define, if exists, the unique y having the largest right overlap with x.
    if y==None: return 0                                    # if no unique right extension, finished, x is not right extensible. 
    y_= kc.get_reverse_sr(y)
    xprime_, dontcare = right_unique_extention(SR,y_)
    if xprime_==None: return 0
    
    # if x != kc.get_reverse_sr(xprime_):        #
        # print(x)
        # print(y)
        # print(len_u)
        # print(kc.get_reverse_sr(xprime_))
        # print(dontcare)
        # print(SR.get_lists_starting_with_given_prefix([46315]))
    assert(x==kc.get_reverse_sr(xprime_))

    # ***** FUSION *****
    SR.remove(x)
    SR.remove(xprime_)
    SR.remove(y)
    SR.remove(kc.get_reverse_sr(y))
    new=x+y[len_u:]
    SR.sorted_add(new)
    SR.sorted_add(kc.get_reverse_sr(new))
    return 1
    
         
def compaction(SR):
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
        witness = fusion(SR,sr)
        if witness == 1: # a fusion was done
            compacted+=1
    sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\n")
    return SR
    
             
def main():
    '''
    Compaction of set of super reads coded as set of ids of unitigs
    '''
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
    SR=compaction(SR)
    sys.stderr.write("  Compaction of simple paths. Done - nb SR="+ str(len(SR))+"\n")
    
    
    sys.stderr.write("  Print maximal super reads - nb MSR="+ str(len(SR))+"\n")
    kc.print_maximal_super_reads(SR)

if __name__ == "__main__":
     main()  
