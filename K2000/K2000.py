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
import K2000_common as kc
def remove_y_strictly_included_in_x(x,SR):
    ''' remove all y strictly included in x'''
    n=len(x)
    # print ("x",x)
    # print ("SUFFIXES")
    for position_suffix in range(0,n):
        u=x[position_suffix:]
        Y=kc.get_SR_starting_with_given_prefix(SR,u)
        # print(u,":",Y)
        if x in Y: Y.remove(x)
        for y in Y:
            # print ("y",y)
            if len(y)+position_suffix < n and colinear(x,[y],[position_suffix]): 
                SR.remove(y)
                SR.remove(kc.get_reverse_sr(y))
    
    # print ("PREFIXES")
    for ending_prefix in range(1,n):
        u=x[:1]
        Y=kc.get_SR_starting_with_given_prefix(SR,u)
        # print(u,":'",Y)
        if x in Y: Y.remove(x)
        for y in Y:
            # print ("y'",y)
            if len(y)< n and colinear(x,[y],[0]): 
                SR.remove(y)
                SR.remove(kc.get_reverse_sr(y))
        
    
    
    

def remove_strict_inclusions(SR):
    ''' remove all super reads stricly included in any other '''
    SR=kc.unique(SR)
    for sr in SR: remove_y_strictly_included_in_x(sr,SR)
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

def fusion (SR,x):
    ''' Main function. For a given super read x, we find y that overlap x, such that : 
     1/ there exists no other y with bigger overlap
     2/ for this y there exist no other x with bigger overlap
     3/ there exists no other x' having an overlap with y and that does not overlap perfectly with x (left colinearity)
     4/ there exists no other y' having an overlap with x and that does not overlap perfectly with y (right colinearity)
    Returns -1 if x cannot be compacted, 0 if it may be compacted latter, 1 if it was compacted
    '''
    # print ("fusion", x)
    n=len(x)
    x_=kc.get_reverse_sr(x)
    for len_u in range(n-1,0,-1): # ***** CONDITION 1/ *****
        u=x[-len_u:]
        # print("u=",u)
        Y=kc.get_SR_starting_with_given_prefix(SR,u)
        assert( x not in Y) 
        # print ("Y=",Y)
        if len(Y)==0: continue  # No y starting with u

        if len(Y)>1: return -1   # More than one unique y starting with u, for instance y and y'. Knowing that y is not included in y' it means necessary that y and y' are not colinear and that x is thus non compactable. 
       
        y=Y[0]
        
        
        
        # ***** CONDITION 2/ *****
        # check if u is maximal for x: we do this by checking x_ (revcomp of x) that start with a suffix of a y_
        y_=kc.get_reverse_sr(y)
        for len_v in range(n,len_u,-1):
            v=y_[-len_v:]
            X_=kc.get_SR_starting_with_given_prefix(SR,v)
            if y_ in X_: X_.remove(y_)
            # print("v=",v," X_=",X_," len(X_)=",len(X_))
            if len(X_)>0: # at least one  other x have a better overlap with y, we stop the compaction for x, but we keep it as a possible later compactable
                return 0
        
        
       
        # ***** CONDITION 3/ *****
        # get all x' such that LCSP(x',y) in [1,len(u)-1]
        # # get all x' ending with a prefix of u ==
        # get all x' starting with a reverse complement of a prefix of u
        X=[]
        ending_positions=[]
        for size_prefix in range(1,len_u+1):
            prefix_u=u[:size_prefix]
            others = kc.get_SR_starting_with_given_prefix(SR,kc.get_reverse_sr(prefix_u))
            if len(others) >0:
                X+=others
                ending_positions+=[len_u-size_prefix for zz in range(len(others))]
        # print ("X=",X,"ending_positions = ",ending_positions)
        if len(ending_positions)>0 and not colinear(x_,X,ending_positions): return -1# ending positions are starting positions of reversed sequences
        
       
        # ***** CONDITION 4/ *****
        # get all y' such that LCSP(x,y') in [1,len(u)-1]
        # get all y' starting with a suffix of u
        Y=[]
        starting_positions=[]
        for starting_suffix_position in range(0,len_u):
            suffix_u=u[starting_suffix_position:]
            others = kc.get_SR_starting_with_given_prefix(SR,suffix_u)
            if len(others) >0:
                Y+=others
                starting_positions+=[starting_suffix_position for zz in range(len(others))]
        # print ("Y=",Y,"ending_positions = ",starting_positions)
        if len(starting_positions)>0 and not colinear(y,Y,starting_positions): return -1
            
        
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
        SR+=[new]
        new=[kc.get_reverse_sr(new)]
        SR+=new
        SR.sort()
        return 1
         
             
def main():
    '''
    Compaction of set of super reads coded as set of ids of unitigs
    '''
    # SR=[[1,3,4],[14],[4,6],[-6,-1,-2],[4,6,7,9],[9,10,11,12],[4,6],[-13, -12, -11]]
    
    SR=kc.generate_SR(sys.argv[1])
    # SR=[[6063, 3129], [6063, 3129, -4346]]
    kc.add_reverse_SR(SR)
    SR.sort()
    SR=remove_strict_inclusions(SR)
    # print("SR",SR)
    dont_try=[]
    while True:
        modified=False
        for sr in SR:
            if kc.contains(sr,dont_try):
                continue
            witness = fusion(SR,sr)
            if witness == 1: # a fusion was done
                modified=True
            elif witness == -1: # no fusion done and sr no compactable 
                dont_try+=[sr]
                # print ("dont try", sr)
            # if witness == 0 : we do nothing, sr, may be compacted latter. 
        if not modified : break
    kc.print_maximal_super_reads(SR)

if __name__ == "__main__":
     main()  
