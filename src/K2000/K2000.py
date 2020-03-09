#!/usr/bin/env python3
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



# def compaction(SR, unitig_lengths,k,min_conflict_overlap): # TO REMOVE
#     ''' Compaction of all sr in SR
#     If a sr was not compacted, i will never be compacted.
#     If it was compacted, maybe it will be re-compacted later on. However, no need to re-run the fusion on this read as
#      - either I could have been compacted on the left and this was done before or this will be done latter or
#      - it will be right extended later: the 'new' (the compacted sr) sr is higher in the lexicographic order than the original sr (as it is longer), thus its position is higher in the SR data structure, thus it will be seen again later.
#     Note that this may be not true in parallel computations.
#     '''
#
#     checked=0
#     compacted=0
#     n = len(SR)
#     for sr in SR.traverse():
#         if checked%100==0:
#             sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\r")
#         checked+=1
#         witness = fusion(SR,sr,unitig_lengths,k,min_conflict_overlap)
#         if witness == 1: # a fusion was done
#             compacted+=1
#     sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\n")
#     return SR



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




# def get_reverse_sr(x):
#     ''' reverse of a super read x. Example is super read x = [4,2,6], reverse(x) = [-6,-2,-4] '''
#     return [-b for b in x][::-1]



def remove_bulles2(SR,unitig_lengths,k,bulles_c):
    checked=0
    n = len(SR)
    for sr in SR.traverse():
        if checked%100==0: sys.stderr.write("      Removing bulles, "+str(checked)+" checked. Size SR "+str(len(SR))+" %.2f"%(100*checked/n)+"%\r")
        checked+=1
        kc.clean_complex_bulles(SR,sr,unitig_lengths,k,bulles_c)
        kc.clean_complex_bulles(SR,kc.get_reverse_sr(sr),unitig_lengths,k,bulles_c)
    sys.stderr.write("Removing bulles, "+str(checked)+" checked. Size SR "+str(len(SR))+" 100%\n")
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
                        help=" NOT IMPLEMENTED WITH ABUNDANCE. Dead end smaller or equal than this value are removed from the path graph.\n Note that with C>0, size of unitigs has to be computable, thus K2000 needs to know the k value and the unitig length. Thus, with C>0, options -k and --unitig_file  are mandatory. [DEFAULT 0]", default=0)

    parser.add_argument("-u", "--unitig_file", type=str,
                        help=" input fasta file containing unitig sequences. Note that unitig id 1 (as indicated in the paths input file) corresponds to the first unitig. This option is mandatory if -c > 0 or -t > 0, else it's useless")

    parser.add_argument("-k", type=int, dest='k',
                        help="kmer size. This option is mandatory if -c > 0 or -t > 0, else it's useless.")

    parser.add_argument("-b", type=int, dest='b',
                        help=" NOT IMPLEMENTED WITH ABUNDANCE. bulle crush", default=0)

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
    csr=kc.generate_SR(input_file)
    sys.stderr.write("  Load super reads. Done - nb SR="+ str(len(csr.SR))+"\n")

    sys.stderr.write("  Add reverse complements \r")
    csr.add_reverses()
    sys.stderr.write("  Add reverse complements. Done - nb SR="+ str(len(csr.SR))+"\n")
    
    for i in 1,2:
        sys.stderr.write("  Compaction of simple paths, min conflict overlap ="+str(min_conflict_overlap)+" \n")
        csr.compaction(unitig_lengths,k,min_conflict_overlap)
        sys.stderr.write("  Compaction of simple paths. Done - nb SR="+ str(len(csr.SR))+"\n")



    if max_tip>0 : # TODO WITH ABUNDANCES
        sys.stderr.write("\n GRAPH CLEANING NOT IMPLEMENTED WITH ABUNDANCES\n\n")
        sys.exit(1)
        sys.stderr.write("\n TIPPING\n\n")
        for x in range(0, 3):
            sys.stderr.write("  Remove tips of size at most "+str(max_tip)+"\n")
            SR=remove_tips(SR,unitig_lengths,k,max_tip)
            sys.stderr.write("  Remove tips. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Remove redundant overlaps\r")
            remove_redundant_overlaps(SR,unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Remove redundant overlaps. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Compaction2 of simple paths \r")
            csr.compaction(unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Compaction2 of simple paths. Done - nb SR="+ str(len(SR))+"\n")

    if(bulles_c>0): # TODO WITH ABUNDANCES
        sys.stderr.write("\n BULLES CRUSHER NOT IMPLEMENTED WITH ABUNDANCES\n\n")
        sys.exit(1)

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
            SR=remove_bulles2(SR,unitig_lengths,k,bulles_c)
            sys.stderr.write("  Remove bulles2. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Remove redundant overlaps\r")
            remove_redundant_overlaps(SR,unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Remove redundant overlaps. Done - nb SR="+ str(len(SR))+"\n")

            sys.stderr.write("  Compaction2 of simple paths \r")
            SR=compaction(SR, unitig_lengths,k,min_conflict_overlap)
            sys.stderr.write("  Compaction2 of simple paths. Done - nb SR="+ str(len(SR))+"\n")


    sys.stderr.write("  Print maximal super reads\n")
    csr.print_maximal_super_reads()

if __name__ == "__main__":
     main()
