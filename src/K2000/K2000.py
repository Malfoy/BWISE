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
import sorted_list
import argparse





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
                    if not kc.is_palindromic(y): SR.remove(kc.get_reverse_sr(y))







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





def right_unique_extention(SR,sr, unitig_lengths,k,min_conflict_overlap):
    ''' return the unique  possible right sr extension with the largest overlap
        return None if no right extensions or if more than one possible non colinear extensions
        The returned extension can be sr itself
    '''

    n=len(sr)
    #  **** Get the largest right overlap with sr ****
    for len_u in range(n-1,0,-1):
        u=sr[-len_u:]
        Y=SR.get_lists_starting_with_given_prefix(u)
        if len(Y)==0: continue              # No y starting with u
        if len(Y)>1: return None,len_u      # More than one unique y starting with u, for instance y and y'. Knowing that y is not included in y' it means necessary that y and y' are not colinear and that x is thus not right extensible.
        y=Y[0]                              # We found the largest y right overlapping sr.

        if min_conflict_overlap>0:
            lenACGT_u = kc.get_len_ACGT(u,unitig_lengths,k)
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
            if min_conflict_overlap>0:
                lenACGT_suffix_u = kc.get_len_ACGT(suffix_u,unitig_lengths,k)           # TODO: optimize this (can be updated from previous loop pass)
                # if lenACGT_u - lenACGT_suffix_u > 500:  break                        # TODO: this value should be a parameter and maybe a ratio.
                # sys.stderr.write("\n          "+str(lenACGT_u)+"  "+str(lenACGT_suffix_u)+"\n")
                #~ if (lenACGT_suffix_u < min_conflict_overlap) : break
                if lenACGT_u - lenACGT_suffix_u > min_conflict_overlap : break
                # if (lenACGT_u / lenACGT_suffix_u) > 10: break                        # TODO: this value should be a parameter
            # END OF THE GREADY PART
            others = SR.get_lists_starting_with_given_prefix(suffix_u)
            if len(others) >0:
                Y+=others
                starting_positions+=[starting_suffix_position for zz in range(len(others))]

        if len(starting_positions)>0 and not kc.colinear(y,Y,starting_positions): return None,len_u     # more than a unique right extention for x.
        return y,len_u                                                                                  # A unique maximal right extention for x (unless gready: a unique largest extention, smaller possible extentions under the gready threahold)

    # not any right extention found.
    return None,None



def fusion (SR,x, unitig_lengths,k,min_conflict_overlap):
    '''Main function. For a given super read x, we find y that overlap x with the highest overlap, such that :
    1/ there exists no other y' right overlapping x that is not collinear with y
    2/ there exists no other x' left overlapping y that is not collinear with x
    Once done, we compact x and y, and this is finished for x.
    '''

    y,len_u=right_unique_extention(SR,x, unitig_lengths,k,min_conflict_overlap)                 # Define, if exists, the unique y != x having the largest right overlap with x.
    if y==None: return 0                                                                        # if no unique right extension, finished, x is not right extensible.
    if y==x: return 0                                                                           # Do not compact x with itself, else, enter an infinite loop
    y_= kc.get_reverse_sr(y)                                                                    # what are left extentions of y, we right extend its reverse complement.
    if y_ == x: return 0                                                                        # Do not compact x with its own reverse complement.
    xprime_, dontcare = right_unique_extention(SR,y_, unitig_lengths,k,min_conflict_overlap)    # Define, if exists, the unique xprime_ (!= y_) having the largest right overlap with y_.
    if xprime_==None: return 0                                                                  # if no unique left extension of the unique right extention of x, finished, x is not right extensible.


    if min_conflict_overlap:                                                                    # if gready its possible* (see below) that y_ is extended with something else than x_
        if x!=kc.get_reverse_sr(xprime_): return 0                                              # In this case, do not make the x/y fusion.
    #~ else:
        #~ assert(x==kc.get_reverse_sr(xprime_))                                                   # if not gready the unique right extension of reverse complement of x is x_

    # * "if gready its possible"* eg:
    # x  -----------
    # y           ----------------
    # x'   /------------------         # x is extended only with y and y_ is extended only with x'_ (as the x/y overlap is too small
    # Note that latter: y will be compacted with x' (when the entry of the fusing function is y'_


    # ***** FUSION *****
    # 1/ remove x and its reverse complement if not palindromic
    # 2/ if y is not x (x==y is possible if x is repeated 2,2,2,2 for instance or if prefix = suffix (1,2,1 for instance)), remove y and its reverse complement if not palindromic
    # 3/ create the new xy SR and add it (sorted fashion)

    # 1
    SR.remove(x)
    if not kc.is_palindromic(x):   SR.remove(xprime_)                                           # Here xprime_ == x_ (we avoid to compute x_ by re-using the xprime_)

    # 2
    if x != y:
        SR.remove(y)
        if not kc.is_palindromic(y):   SR.remove(kc.get_reverse_sr(y))

    # 3
    new=x+y[len_u:]
    SR.sorted_add(new)
    if not kc.is_palindromic(new): SR.sorted_add(kc.get_reverse_sr(new))

    # we made a compaction, return 1.
    return 1



def remove_redundant_overlap(SR,sr, unitig_lengths,k,min_conflict_overlap):
    right,len_u_right   =right_unique_extention(SR,sr,                      unitig_lengths,k,min_conflict_overlap)
    if right==None: return
    left,len_u_left     =right_unique_extention(SR,kc.get_reverse_sr(sr),   unitig_lengths,k,min_conflict_overlap)
    if left==None:  return
    if len_u_left+len_u_right > len(sr):                # The overlap is larger than sr. sr is redundant.
        SR.remove(sr)
        if not kc.is_palindromic(sr): SR.remove(kc.get_reverse_sr(sr))



def remove_redundant_overlaps(SR,unitig_lengths,k,min_conflict_overlap):
    for sr in SR.traverse():
        remove_redundant_overlap(SR,sr, unitig_lengths,k,min_conflict_overlap)



def compaction(SR, unitig_lengths,k,min_conflict_overlap):
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
        witness = fusion(SR,sr,unitig_lengths,k,min_conflict_overlap)
        if witness == 1: # a fusion was done
            compacted+=1
    sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\n")
    return SR



def remove_tips(SR,unitig_lengths,k,maxTip):
    checked=0
    n = len(SR)
    for sr in SR.traverse():
        if checked%100==0: sys.stderr.write("      Removing tips, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\r")
        checked+=1
        if kc.get_len_ACGT(sr,unitig_lengths,k) < maxTip and kc.to_clean(SR,sr):
            SR.remove(sr)
            if not kc.is_palindromic(sr): SR.remove(kc.get_reverse_sr(sr))
        #~ if kc.get_len_ACGT(sr,unitig_lengths,k) < maxTip/10 and kc.is_a_dead_end(SR,sr):
            #~ SR.remove(sr)
            #~ if not kc.is_palindromic(sr): SR.remove(kc.get_reverse_sr(sr))
    sys.stderr.write("      Removing tips, "+str(checked)+" checked. Size SR "+str(len(SR))+" 100%\n")
    return SR


def remove_bulles(SR):
    checked=0
    n = len(SR)
    for sr in SR.traverse():
        if checked%100==0: sys.stderr.write("      Removing bulles, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\r")
        checked+=1
        kc.clean_parallel_contigs(SR,sr)

    sys.stderr.write("      Removing bulles, "+str(checked)+" checked. Size SR "+str(len(SR))+" 100%\n")
    return SR


def remove_bulles2(SR):
    checked=0
    n = len(SR)
    for sr in SR.traverse():
        if checked%100==0: sys.stderr.write("      Removing bulles, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\r")
        checked+=1
        kc.clean_complex_bulles(SR,sr)
        kc.clean_complex_bulles(SR,get_reverse_sr(sr))
    sys.stderr.write("      Removing bulles, "+str(checked)+" checked. Size SR "+str(len(SR))+" 100%\n")
    return SR



def remove_dust(SR,unitig_lengths,k,maxDust):
    for sr in SR.traverse():
        if kc.get_len_ACGT(sr,unitig_lengths,k) < maxDust:
            SR.remove(sr)
            if not kc.is_palindromic(sr): SR.remove(kc.get_reverse_sr(sr))
    return SR


def main():
    '''
    Compaction of set of super reads coded as set of ids of unitigs
    '''

    parser = argparse.ArgumentParser(description='Compaction of set of super reads coded as set of ids of unitigs.')
    parser.add_argument("input_file", type=str,
                        help="input file containing dbg paths as a list of unitig ids, eg. on line looks like \"-1;24;198;\"" )

    parser.add_argument("-c", "--min_conflict_overlap", type=int, dest='c',
                        help="Minimal conflict overlap. \n\t With c=0: K2000 is exact. \n\t With c>0: K2000 becomes greedy, in this case if a path A could be extended either by B or C and B and C are not collinear, then if the size of the overlap (A,B)>c and the size of the overlap (A,C)<c, then compact A-B but not A-C. If both overlaps are bigger than c or both smaller than c, no compaction is made. \n Note that with c>0, size of unitigs has to be computable, thus K2000 needs to know the k value and the unitig length. Thus, with c>0, options -k and --unitig_file  are mandatory. [DEFAULT 0]", default=0)

    parser.add_argument("-t", "--max_tip", type=int, dest='t',
                        help=" Dead end smaller or equal than this value are removed from the path graph.\n Note that with C>0, size of unitigs has to be computable, thus K2000 needs to know the k value and the unitig length. Thus, with C>0, options -k and --unitig_file  are mandatory. [DEFAULT 0]", default=0)

    parser.add_argument("-u", "--unitig_file", type=str,
                        help=" input fasta file containing unitig sequences. Note that unitig id 1 (as indicated in the paths input file) corresponds to the first unitig. This option is mandatory if -c > 0 or -t > 0, else it's useless")

    parser.add_argument("-k", type=int, dest='k',
                        help="kmer size. This option is mandatory if -c > 0 or -t > 0, else it's useless.")

    parser.add_argument("-b", type=int, dest='b',
                        help="bulle crush", default=0)

    args = parser.parse_args()
    input_file=str(args.input_file)
    min_conflict_overlap=args.c
    max_tip=args.t
    k=args.k
    bulles_c=args.b
    unitig_file=args.unitig_file


    sys.stderr.write("** This is K2000. Option reviews     **\n")
    sys.stderr.write("\t input_file: "+            input_file+                  "\n")
    sys.stderr.write("\t min_conflict_overlap: "+  str(min_conflict_overlap)+   "\n")
    sys.stderr.write("\t max tips: "+              str(max_tip)+                "\n")
    sys.stderr.write("\t k: "+                     str(k)+                      "\n")
    sys.stderr.write("\t unitig_file: "+           str(unitig_file)+                 "\n")
    sys.stderr.write("** This is K2000. Computation starts **\n")

    # User needs to remove tips or uses a gready approach. In this case we need to check that k and unitig lengths are correctly informed, then we load unitigs lengths.
    unitig_lengths=None
    if max_tip>0 or min_conflict_overlap>0:
        if k==None or unitig_file==None:
            sys.stderr.write("ERROR: k option and unitig_file must be informed if using max_tip>0 or min_conflict_overlap>0. Exit\n")
            sys.exit(1)
        sys.stderr.write("  Load unitig lengths \n")
        unitig_lengths = kc.load_unitig_lengths (unitig_file)


    sys.stderr.write("  Load super reads \r")
    SR=kc.generate_SR(input_file)
    sys.stderr.write("  Load super reads. Done - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Add reverse complements \r")
    kc.add_reverse_SR(SR)
    sys.stderr.write("  Add reverse complements. Done - nb SR="+ str(len(SR))+"\n")

    #~ sys.stderr.write("  Conserve only unique super reads \r")
    #~ SR.unique()
    #~ sys.stderr.write("  Conserve only unique super reads. Done - nb SR="+ str(len(SR))+"\n")

    #~ sys.stderr.write("  Sorting \r")
    #~ SR.sort()
    #~ sys.stderr.write("  Sorting. Done - nb SR="+ str(len(SR))+"\n")

    #~ sys.stderr.write("  Remove strict inclusions\r")
    #~ SR=remove_strict_inclusions(SR)
    #~ sys.stderr.write("  Remove strict inclusions. Done - nb SR="+ str(len(SR))+"\n")

    #~ if max_tip>0:
        #~ sys.stderr.write("  Remove tips of size at most "+str(max_tip)+"\n")
        #~ SR=remove_tips(SR,unitig_lengths,k,max_tip)
        #~ sys.stderr.write("  Remove tips. Done - nb SR="+ str(len(SR))+"\n")


    sys.stderr.write("  Compaction of simple paths, min conflict overlap ="+str(min_conflict_overlap)+" \n")
    SR=compaction(SR, unitig_lengths,k,min_conflict_overlap)
    sys.stderr.write("  Compaction of simple paths. Done - nb SR="+ str(len(SR))+"\n")

    sys.stderr.write("  Compaction of simple paths, min conflict overlap ="+str(min_conflict_overlap)+" \n")
    SR=compaction(SR, unitig_lengths,k,min_conflict_overlap)
    sys.stderr.write("  Compaction of simple paths. Done - nb SR="+ str(len(SR))+"\n")


    if max_tip>0 :
        sys.stderr.write("\n TIPPING\n\n")
        for x in range(0, 3):
            sys.stderr.write("  Remove tips of size at most "+str(max_tip)+"\n")
            SR=remove_tips(SR,unitig_lengths,k,max_tip)
            sys.stderr.write("  Remove tips. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Remove redundant overlaps\r")
            remove_redundant_overlaps(SR,unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Remove redundant overlaps. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Compaction2 of simple paths \r")
            SR=compaction(SR, unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Compaction2 of simple paths. Done - nb SR="+ str(len(SR))+"\n")

    if(bulles_c>0):
    #~ if(True):
        sys.stderr.write("\n BULLES CRUSH\n\n")
        for x in range(0, 3):
            sys.stderr.write("  Remove bulles \n")
            SR=remove_bulles(SR)
            sys.stderr.write("  Remove bulles. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Remove redundant overlaps\r")
            remove_redundant_overlaps(SR,unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Remove redundant overlaps. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Compaction2 of simple paths \r")
            SR=compaction(SR, unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Compaction2 of simple paths. Done - nb SR="+ str(len(SR))+"\n")


            sys.stderr.write("  Remove bulles2 \n")
            SR=remove_bulles2(SR)
            sys.stderr.write("  Remove bulles2. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Remove redundant overlaps\r")
            remove_redundant_overlaps(SR,unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Remove redundant overlaps. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Compaction2 of simple paths \r")
            SR=compaction(SR, unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Compaction2 of simple paths. Done - nb SR="+ str(len(SR))+"\n")


    sys.stderr.write("  Print maximal super reads\n")
    kc.print_maximal_super_reads(SR)

if __name__ == "__main__":
     main()
