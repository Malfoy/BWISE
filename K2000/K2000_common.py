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
    
def compare (tuple1,tuple2):
    '''tuple2 always smaller than tuple1, no verification made'''
    '''
    if tuple1 starts with tuple2: return 0
    if tuple1<tuple2: return -1
    if tuple1>tuple2: return 1
    '''
    tmp_tuple1=tuple1[0:len(tuple2)]
    if tmp_tuple1 == tuple2: return  0
    if tmp_tuple1 <  tuple2: return -1
    return 1
    

def get_SR_starting_with_given_prefix(SR, prefix):
    ''' given a prefix of a super read, return all super reads in SR starting with this prefix
    Dichotomy. log(|SR|) comparisons
    O(|prefix|*log(|SR|))
    '''
    
    start=0
    n=len(SR)
    stop=n
    
   
    
    while start <=stop: 
        middle = start+int((stop-start)/2)
        if middle>=n: break
        tuple1 = SR[middle]
        cmp_val = compare(tuple1,prefix)
        if cmp_val == -1: # prefix may start in the remaining part of the SR array
            start = middle+1
        if cmp_val ==  1: # prefix may start in the starting part of the SR array
            stop = middle-1
        if cmp_val == 0: # we found a tuple starting with the prefix. We need to check other tuples starting with prefix before and after in the array. 
            res=[tuple1]
            i=middle-1
            while i>0 and compare(SR[i],prefix)==0:
                res.append(SR[i])
                i-=1
            i=middle+1
            while i<n and compare(SR[i],prefix)==0:
                res.append(SR[i])
                i+=1
            return res
    return []
            


def contains (sr,SR):
    ''' if sr is in SR: return True, else return False.
    faster (O(log(|SR|)) that the 'in' function from non ordered array (O(|SR|))
    '''
    X=get_SR_starting_with_given_prefix(SR,sr)
    if sr in X: return True # O(|X|) which can be considered as constant
    return False

def generate_SR(file_name):
    ''' Given an input file storing super reads, store them in the SR array'''
    unitig_file = open(file_name, 'r')
    SR=[]
    for line in unitig_file:
        
        line = line.rstrip()[:-1].split(';')
        sr=[]
        for unitig_id in line: 
            
            sr=sr+[int(unitig_id)]
        SR+=[sr]
    return SR
    
def add_reverse_SR(SR):
    ''' For all super reads in SR, we add there reverse in SR 
    This double the SR size.
    We don't check that this does not create any duplicates'''
    SR_=[]
    for sr in SR:
        SR_+=[get_reverse_sr(sr)]
    for sr in SR_:
        SR+=[sr]
        

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
    for sr in SR:
     if is_canonical(sr):
         if len(sr)==1:
             print (str(sr[0])+";")
         else:
             for x in sr: print (str(x)+";", end="")
             print ()
             



                
def unique(s): # from http://code.activestate.com/recipes/52560-remove-duplicates-from-a-sequence/
    """Return a list of the elements in s, but without duplicates.

    For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
    unique("abcabc") some permutation of ["a", "b", "c"], and
    unique(([1, 2], [2, 3], [1, 2])) some permutation of
    [[2, 3], [1, 2]].

    For best speed, all sequence elements should be hashable.  Then
    unique() will usually work in linear time.

    If not possible, the sequence elements should enjoy a total
    ordering, and if list(s).sort() doesn't raise TypeError it's
    assumed that they do enjoy a total ordering.  Then unique() will
    usually work in O(N*log2(N)) time.

    If that's not possible either, the sequence elements must support
    equality-testing.  Then unique() will usually work in quadratic
    time.
    """

    n = len(s)
    if n == 0:
        return []

    # Try using a dict first, as that's the fastest and will usually
    # work.  If it doesn't work, it will usually fail quickly, so it
    # usually doesn't cost much to *try* it.  It requires that all the
    # sequence elements be hashable, and support equality comparison.
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()

    # We can't hash all the elements.  Second fastest is to sort,
    # which brings the equal elements together; then duplicates are
    # easy to weed out in a single pass.
    # NOTE:  Python's list.sort() was designed to be efficient in the
    # presence of many duplicate elements.  This isn't true of all
    # sort functions in all languages or libraries, so this approach
    # is more effective in Python than it may be elsewhere.
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]

    # Brute force is all that's left.
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u
