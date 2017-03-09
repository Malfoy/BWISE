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

    # pool_size = 8  # your "parallelness"
    # pool = Pool(pool_size)

    for sr in SR.traverse():
        if checked%100==0: sys.stderr.write("Removing inclusions, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\r")
        checked+=1
        # pool.apply_async(remove_y_strictly_included_in_x(sr,SR,to_remove))

        remove_y_strictly_included_in_x(sr,SR)#,to_remove)

    sys.stderr.write("Removing inclusions, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\n")
    # SR_2=[]
  #   for i in range(len(to_remove)):
  #       if not to_remove[i]:
  #           SR_2.append(SR[i])
  #   del SR
  #   return SR_2
    return SR
    
    


    
def colinear(x,X,starting_positions):
    ''' Check that all sr in X are colinear with x
    For each sr in X, one knows its starting position on x, with table starting_positions'''
    # print ("colinear", x,X,starting_positions)
    for i in range(len(X)):
        other = X[i]
        # print ("other",other)
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
        # print ("u is", u)
        Y=SR.get_lists_starting_with_given_prefix(u)
        # print ("Y is", Y)
        if sr in Y: Y.remove(sr)    # possible if x is repeated 2,2,2,2 for instance.
        if len(Y)==0: continue      # No y starting with u
        if len(Y)>1: return None,len_u    # More than one unique y starting with u, for instance y and y'. Knowing that y is not included in y' it means necessary that y and y' are not colinear and that x is thus right extensible. 
        y=Y[0]                      # We found the largest y right overlapping sr.
    
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
    Once done, we compact x and y, and 
    '''
    
    y,len_u=right_unique_extention(SR,x)                    # Define, if exists, the unique y having the largest right overlap with x.
    if y==None: return 0                                    # if no unique right extension, finished, x is not right extensible. 
    y_= kc.get_reverse_sr(y)
    xprime_, dontcare = right_unique_extention(SR,y_)
    if xprime_==None: return 0
    
    if x != kc.get_reverse_sr(xprime_):        #
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
    
    
    
        
    
    
def fusion_old (SR,x):
    ''' Main function. For a given super read x, we find y that overlap x, such that : 
     1/ there exists no other y with bigger overlap
     2/ there exists no other x' having an overlap with y and that does not overlap perfectly with x (left colinearity)
     3/ there exists no other y' having an overlap with x and that does not overlap perfectly with y (right colinearity)
    Returns -1 if x cannot be compacted, 0 if it may be compacted latter, 1 if it was compacted
    '''
    # print ("fusion", x)
    n=len(x)
    # print("x is ", x)
    x_=kc.get_reverse_sr(x)
    for len_u in range(n-1,0,-1): # ***** CONDITION 1/ *****
        # print ("len(u)", len_u)
        u=x[-len_u:]
        # print("u=",u)
        Y=SR.get_lists_starting_with_given_prefix(u)
        if x in Y: Y.remove(x) # possible if x is repeated 2,2,2,2 for instance.
        # print ("Y=",Y)
        if len(Y)==0: continue  # No y starting with u

        if len(Y)>1: return 0   # More than one unique y starting with u, for instance y and y'. Knowing that y is not included in y' it means necessary that y and y' are not colinear and that x is thus non compactable. 
       
        y=Y[0]
        

        
        
       
        # ***** CONDITION 2/ *****
        # get all x' such that LCSP(x',y) in [1,len(u)-1]
        # # get all x' ending with a prefix of u ==
        # get all x' starting with a reverse complement of a prefix of u
        X=[]
        ending_positions=[]
        for size_prefix in range(1,len_u+1):
            prefix_u=u[:size_prefix]
            others = SR.get_lists_starting_with_given_prefix(kc.get_reverse_sr(prefix_u))
            if len(others) >0:
                X+=others
                ending_positions+=[len_u-size_prefix for zz in range(len(others))]
        # print ("X=",X,"ending_positions = ",ending_positions)
        if len(ending_positions)>0 and not colinear(x_,X,ending_positions): return 0# ending positions are starting positions of reversed sequences
        
       
        # ***** CONDITION 3/ *****
        # get all y' such that LCSP(x,y') in [1,len(u)-1]
        # get all y' starting with a suffix of u
        Y=[]
        starting_positions=[]
        for starting_suffix_position in range(0,len_u):
            suffix_u=u[starting_suffix_position:]
            others = SR.get_lists_starting_with_given_prefix(suffix_u)
            if len(others) >0:
                Y+=others
                starting_positions+=[starting_suffix_position for zz in range(len(others))]
        # print ("Y=",Y,"ending_positions = ",starting_positions)
        if len(starting_positions)>0 and not colinear(y,Y,starting_positions): return 0
            
        
        # print ("fusionning ",x,"with",y,", generating", x+y[len_u:])
        # ***** FUSION *****
        # print ("remove",x)
        SR.remove(x)
        # print ("remove",x_)
        SR.remove(x_)
        # print ("remove",y)
        SR.remove(y)
        SR.remove(kc.get_reverse_sr(y))
        # print ("remove",kc.get_reverse_sr(y))
        new=x+y[len_u:]
        SR.sorted_add(new)
        SR.sorted_add(kc.get_reverse_sr(new))
        return 1
         

             
def main():
    '''
    Compaction of set of super reads coded as set of ids of unitigs
    '''
    # SR=[[1,3,4],[14],[4,6],[-6,-1,-2],[4,6,7,9],[9,10,11,12],[4,6],[-13, -12, -11]]

    sys.stderr.write("Load super reads \n")
    SR=kc.generate_SR(sys.argv[1])
    # SR.unique()
    sys.stderr.write("Add reverse complements "+ str(len(SR))+"\n")
    kc.add_reverse_SR(SR)
    SR.unique()
    sys.stderr.write("Sorting"+ str(len(SR))+"\n")
    SR.sort()
    sys.stderr.write("Remove strict inclusions "+ str(len(SR))+"\n")
    SR=remove_strict_inclusions(SR)
    
    # print("SR",SR)
    # dont_try = sorted_list.sorted_list()
    checked=0
    compacted=0
    n = len(SR)
    for sr in SR.traverse():
        if checked%100==0: 
            sys.stderr.write("Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" nodes were compacted\r")
        checked+=1
        # if dont_try.contains(sr):
            # continue
        witness = fusion(SR,sr)
        if witness == 1: # a fusion was done
            compacted+=1
        # elif witness == -1: # no fusion done and sr no compactable
            # dont_try.sorted_add(sr)
            # print ("dont try", sr)
        
        # if witness == 0 : we do nothing, sr, may be compacted latter. 
    sys.stderr.write("Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" nodes were compacted\n")
    
    sys.stderr.write("Print results "+ str(len(SR))+"\n")
    kc.print_maximal_super_reads(SR)

if __name__ == "__main__":
     main()  
