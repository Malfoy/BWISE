#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
'''
Compaction of super reads. Non ambiguous paths from a set of super reads are compacted
Resulting super reads are called MSR for Maximal Super Reads
Common file
@author (except for the 'unique' function) pierre peterlongo pierre.peterlongo@inria.fr
'''            

import sys
import getopt
import sorted_list
    
    
    
# Create structures used by dict and array methods.
basecomplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
def complement(s): 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)
    

def reverse_complement(seq):
    return ''.join(basecomplement[c] for c in seq[::-1])
    
def get_reverse_sr(x):
    ''' reverse of a super read x. Example is super read x = [4,2,6], reverse(x) = [-6,-2,-4] '''
    return [-b for b in x][::-1]
    
# def compare (tuple1,tuple2):
#     '''
#     if tuple1 starts with tuple2: return 0
#     if tuple1<tuple2: return -1
#     if tuple1>tuple2: return 1
#     '''
#
#     tmp_tuple1=tuple1[0:len(tuple2)]
#     if tmp_tuple1 < tuple2: return    -1
#     if tmp_tuple1 > tuple2: return     1
#     return                             0
#
#
# def get_SR_starting_with_given_prefix(SR, prefix):
#     ''' given a prefix of a super read, return all super reads in SR starting with this prefix
#     Dichotomy. log(|SR|) comparisons
#     O(|prefix|*log(|SR|))
#     '''
#
#     start=0
#     n=len(SR)
#     stop=n
#
#
#
#     while start <=stop:
#         middle = start+int((stop-start)/2)
#         if middle>=n: break
#         tuple1 = SR[middle]
#         cmp_val = compare(tuple1,prefix)
#         if cmp_val == -1:   # prefix may start in the remaining part of the SR array
#             start = middle+1
#             continue
#         if cmp_val ==  1:   # prefix may start in the starting part of the SR array
#             stop = middle-1
#             continue
#         # if cmp_val == 0:    # we found a tuple starting with the prefix. We need to check other tuples starting with prefix before and after in the array.
#         res=[tuple1]
#         i=middle-1
#         while i>-1 and compare(SR[i],prefix)==0:
#             res.append(SR[i])
#             i-=1
#         i=middle+1
#         while i<n and compare(SR[i],prefix)==0:
#             res.append(SR[i])
#             i+=1
#         return res
#     return []
            

# def contains (sr,SR):
#     ''' if sr is in SR: return True, else return False.
#     faster (O(log(|SR|)) that the 'in' function from non ordered array (O(|SR|))
#     '''
#     X=get_SR_starting_with_given_prefix(SR,sr)
#     if sr in X: return True # O(|X|) which can be considered as constant
#     return False

def generate_SR(file_name):
    ''' Given an input file storing super reads, store them in the SR array'''
    sr_file = open(file_name, 'r')
    sl = sorted_list.sorted_list()
    for line in sr_file:
        
        line = line.rstrip()[:-1].split(';')
        sr=[]
        for unitig_id in line: 
            sr_val=int(unitig_id)
            sr=sr+[sr_val]    
        sl.add(sr)
    return sl
    
def add_reverse_SR(SR):
    ''' For all super reads in SR, we add there reverse in SR 
    This double the SR size.
    We don't check that this does not create any duplicates'''
    SR_ = sorted_list.sorted_list()
    for sr in SR.traverse():
        SR_.add(get_reverse_sr(sr))
    for sr in SR_.traverse():
        SR.add(sr)
    return SR
        

def canonical(sr):
    ''' return the canonical representation of a sr (the smallest version of a sr and its reverse complement'''

    if len(sr)==1:
        if sr[0]>0:
            return [sr]
        else:
            return [-sr]
    sr_=get_reverse_sr(sr)
    if sr>sr_:  return sr
    else :      return sr_
    
def is_canonical(sr):

    ''' return True if the canonical representation of sr is itself'''
    if len(sr)==1:
        if sr[0]>0: return True
        else:       return False
    if sr>get_reverse_sr(sr):    return True
    else:                           return False
    
         
def print_maximal_super_reads(SR):
    '''print all maximal super reads as a flat format'''
    for sr in SR.traverse():
     if is_canonical(sr):
         if len(sr)==1:
             print (str(sr[0])+";")
         else:
             for unitig_id in sr: 

                 # if unitig_id>=0: unitig_id-=1
   #               else: unitig_id+=1
                 print (str(unitig_id)+";", end="")
             print ()
             

            
def load_unitigs(file_name):
    unitigs=[]
    file=open(file_name,'r')
    for line in file:
        line=line.rstrip()
        if line[0]=='>': continue
        unitigs+=[line]
    return unitigs


# SR=generate_SR("test.txt")
# SR.unique()
# SR=add_reverse_SR(SR)
# SR.unique()
# for sr in SR.traverse():
#     print (sr)
