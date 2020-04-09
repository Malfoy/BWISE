#!/usr/bin/env python3
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


class counted_super_reads():
    def __init__(self):
        self.SR=sorted_list.sorted_list()
        self.abundances = {}
    

    
    def add(self,sr,abundance):
        t_sr = tuple(sr)
        if t_sr not in self.abundances:
            self.SR.add(sr)
            self.abundances[t_sr]=abundance
        else:
            self.abundances[t_sr]+=abundance
            
    def sort(self):
        self.SR.sort()
        
    def index_nodes(self):
        self.SR.index_nodes()
        
    def show_right_edges (self,x,id_x,unitigs,k):
        ''' Main function. For a given super read x, we find y that overlap x, and we print the links in a GFA style:
        L    11    +    12    -    overlap size
        Note that one treat x only if its canonical.
        Four cases :
        1/ x overlaps y, with y canonical. One prints x + y + blabla
        2/ x overlaps y_, with y_ non canonical. In this case, y_ overlaps x. One of the two solutions has to be chosen. We chose min(idx,idy) (with idx,idy being the ids of the MSR x,y in SR) One searches the id of y, and one prints x + y - blabla.
        3/ x_ overlaps y. same as 2.
        4/ x_ overlaps y_. We do nothing, this case is treated when the entry of the function is y that thus overlaps x.
        WARNING: here x and each msr in MSR contain as last value its unique id.
        '''
        x=x[:-1]                                # remove the x_id from the x msr
        if not is_canonical(x): return
        n=len(x)

        # CASES 1 AND 2
        strandx='+'
        # print ("x is", x)
        for len_u in range(1,n): # for each possible x suffix
            u=x[-len_u:]
            # print ("u is", u)
            Y=self.SR.get_lists_starting_with_given_prefix(u)
            # if x in Y: Y.remove(x)            # we remove x itself from the list of y : note that it should not occur.
            # print (x)
            # assert(x not in Y)
            if len(Y)==0: continue              # No y starting with u
            for y in Y:
                # detect the y strand
                # CASE 1/
                if is_canonical(y[:-1]):     # remove the last value that corresponds to the node id
                    strandy ='+'
                    # id_y=indexed_nodes.index(y)                               # get the id of the target node
                    id_y=get_msr_id(y)                                           # last value is the node id, here the id of the target node
                # CASE 2/
                else:
                    strandy='-'
    #                id_y = indexed_nodes.index(kc.get_reverse_msr(y))
                    id_y=get_reverse_msr_id(y,self.SR)                                # find the reverse of list y in self.SR to grab its id.
                    if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y.
                size_super_read = get_size_super_read_in_u(u,unitigs,k)
                # print the edges
                print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(size_super_read)+"M")


        # CASES 3 AND 4
        strandx='-'
        x_=get_reverse_sr(x)
        for len_u in range(1,n): # for each possible x suffix
            u=x_[-len_u:]
            Y=self.SR.get_lists_starting_with_given_prefix(u)
            # assert(x_ not in Y)
            if len(Y)==0: continue  # No y starting with u
            for y in Y:
                if is_canonical(y[:-1]): # CASE 3
                    strandy ='+'
                   # id_y=indexed_nodes.index(y)                                # get the id of the target node
                    id_y=get_msr_id(y)                                           # last value is the node id, here the id of the target node
                    # we determine min(id_x,id_y)
                    if id_x>id_y: continue # x_.y is the same as y_.x. Thus we chose one of them. By convention, we print x_.y if x<y.
                    size_super_read = get_size_super_read_in_u(u,unitigs,k)
                    print ("L\t"+str(id_x)+"\t"+strandx+"\t"+str(id_y)+"\t"+strandy+"\t"+str(size_super_read)+"M") # note that strand x is always '-' and strandy is always '+' in this case.

    #            else: continue # CASE 4, nothing to do.
    
    
    def print_GFA(self, unitigs, size_overlap):
        '''print canonical unitigs
        WARNING: here each msr in MSR contains as last value its unique id.
        '''
        nb_errors=0
        ### PRINT NODES
        sys.stderr.write("Print GFA Nodes\n")
        for indexed_msr in self.SR.traverse():
            node_id = get_msr_id(indexed_msr)                                 # last value is the node id
            msr = indexed_msr[:-1]                                            # remove the last value that corresponds to the node id
            if node_id%100==0: sys.stderr.write("\t%.2f"%(100*node_id/len(self.SR))+"%\r")
            if not is_canonical(msr): continue

            print ("S\t"+str(node_id)+"\t", end="")
            print_first_kmer=True
            previous_id=-1
            previous_overlap=""                                         #used only to assert a good k-1 overlap.
            size_msrACGT=0
            msrACGT=""
            for unitig_id in msr:
                reverse=False
                # print ("\n",str(unitig_id-1))
                if unitig_id<0:                                         #the sequence is indicated as reverse complemented. Note that unitig ids start with 1, thus the -0 problem does not appear.
                    reverse=True
                    unitig_id=-unitig_id
                unitig=unitigs[unitig_id-1]                             # grab the good unitig. Ids starts at 1 (avoiding the -0=0 problem). Thus unitig i corresponds to unitigs[i-1]
                # sum_coverage+=(unitig_coverages[unitig_id]*len(unitig)) # The unitig had been seen unitig_coverages[unitig_id] times. As downstreat visualization tool as bandage divides by the sequence size, we multiply the coverage to the unitig by the size of the unitig
                if reverse: unitig=reverse_complement(unitig)        #reverse the untig if necessary
                if previous_overlap != "":                              # overlap validation
                    if(unitig[:size_overlap] != previous_overlap):      # overlap validation
                        nb_errors+=1
                        # sys.stderr.write("\n WARNING unitigs ids "+ str(previous_id)+" and "+str(unitig_id)+" do not overlap, while they do in "+str(msr)+"\n")
                        # sys.stderr.write("Suffix "+previous_overlap+" size="+str(len(previous_overlap))+"\n")
                        # sys.stderr.write("Prefix "+unitig[:size_overlap]+"\n")
                previous_overlap = unitig[-size_overlap:]               # store the suffix of size k-1 to check the next overlap
                previous_id = unitig_id
                if not print_first_kmer: unitig=unitig[size_overlap:]   # remove the k-1 overlap
                print_first_kmer=False                                  #do not print first kmer for next unitigs
                # print (unitig,end="")
                msrACGT+=unitig
            # size_msrACGT=len(msrACGT)
            # coverage=number_mapped_sr[node_id]/float(size_msrACGT)
            len_msrACGT = len(msrACGT)
            fieldRC = len_msrACGT*self.get_abundance(msr)
            print (msrACGT+"\tLN:i:"+str(len_msrACGT)+"\tRC:i:"+str(fieldRC))#+"\t#AVG:f:%.2f"%(coverage))



        sys.stderr.write("\t100.00% -- "+ str(nb_errors)+" error(s)\n" )
        
        ### EDGES
        sys.stderr.write("Print GFA Edges\n")
        for msr in self.SR.traverse():
            x_id = get_msr_id(msr)                                         # last value is the node id
            if x_id%100==0: sys.stderr.write("\t%.2f"%(100*x_id/len(self.SR))+"%\r")
            self.show_right_edges(msr,x_id,unitigs,size_overlap+1)
        sys.stderr.write("\t100.00%\n")
    
    def get_abundance(self,sr):
        ''' return the abundance of a sr 
        the abundance is stored either as the sr or as its reverse complement
        '''
        t_sr = tuple(sr)
        if t_sr in self.abundances: return self.abundances[t_sr]
        t_rsr = tuple(get_reverse_sr(sr))
        assert t_rsr in self.abundances
        return self.abundances[t_rsr]
            
    def add_reverses(self):
        ''' For all super reads in SR, we add there reverse in SR
        This double the SR size, unless there are palindromes ([1,-1] for instance). Those are not added.
        We don't check that this does not create any duplicates'''
        SR_ = sorted_list.sorted_list()
        for sr in self.SR.traverse():
            if not is_palindromic(sr):
                rsr = get_reverse_sr(sr)
                SR_.add(rsr)
        for sr in SR_.traverse():
            self.SR.add(sr)
            
    def right_unique_extention(self, sr, unitig_lengths,k,min_conflict_overlap):
        ''' return the unique  possible right sr extension with the largest overlap
            return None if no right extensions or if more than one possible non colinear extensions
            The returned extension can be sr itself
        '''
        n=len(sr)
        #  **** Get the largest right overlap with sr ****
        for len_u in range(n-1,0,-1):
            u=sr[-len_u:]
            Y=self.SR.get_lists_starting_with_given_prefix(u)
            if len(Y)==0: continue              # No y starting with u
            if len(Y)>1: return None,len_u      # More than one unique y starting with u, for instance y and y'. Knowing that y is not included in y' it means necessary that y and y' are not colinear and that x is thus not right extensible.
            y=Y[0]                              # We found the largest y right overlapping sr.

            if min_conflict_overlap>0:
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
                if min_conflict_overlap>0:
                    lenACGT_suffix_u = get_len_ACGT(suffix_u,unitig_lengths,k)           # TODO: optimize this (can be updated from previous loop pass)
                    # if lenACGT_u - lenACGT_suffix_u > 500:  break                        # TODO: this value should be a parameter and maybe a ratio.
                    # sys.stderr.write("\n          "+str(lenACGT_u)+"  "+str(lenACGT_suffix_u)+"\n")
                    #~ if (lenACGT_suffix_u < min_conflict_overlap) : break
                    if lenACGT_u - lenACGT_suffix_u > min_conflict_overlap : break
                    # if (lenACGT_u / lenACGT_suffix_u) > 10: break                        # TODO: this value should be a parameter
                # END OF THE GREADY PART
                others = self.SR.get_lists_starting_with_given_prefix(suffix_u)
                if len(others) >0:
                    Y+=others
                    starting_positions+=[starting_suffix_position for zz in range(len(others))]

            if len(starting_positions)>0 and not colinear(y,Y,starting_positions): return None,len_u     # more than a unique right extention for x.
            return y,len_u                                                                                  # A unique maximal right extention for x (unless gready: a unique largest extention, smaller possible extentions under the gready threahold)

        # not any right extention found.
        return None,None
    
    def fusion (self, x,unitig_lengths,k,min_conflict_overlap):
        '''Main function. For a given super read x, we find y that overlap x with the highest overlap, such that :
        1/ there exists no other y' right overlapping x that is not collinear with y
        2/ there exists no other x' left overlapping y that is not collinear with x
        Once done, we compact x and y, and this is finished for x.
        '''

        y,len_u=self.right_unique_extention(x, unitig_lengths,k,min_conflict_overlap)                 # Define, if exists, the unique y != x having the largest right overlap with x.
        if y==None: return 0                                                                        # if no unique right extension, finished, x is not right extensible.
        if y==x: return 0                                                                            # Do not compact x with itself, else, enter an infinite loop
        y_= get_reverse_sr(y)                                                                    # what are left extentions of y, we right extend its reverse complement.
        if y_ == x: return 0                                                                        # Do not compact x with its own reverse complement.
        xprime_, dontcare = self.right_unique_extention(y_, unitig_lengths,k,min_conflict_overlap)    # Define, if exists, the unique xprime_ (!= y_) having the largest right overlap with y_.
        if xprime_==None: return 0                                                                    # if no unique left extension of the unique right extention of x, finished, x is not right extensible.


        if min_conflict_overlap:                                                                    # if gready its possible* (see below) that y_ is extended with something else than x_
            if x!=get_reverse_sr(xprime_): return 0                                              # In this case, do not make the x/y fusion.
        #~ else:
            #~ assert(x==.get_reverse_sr(xprime_))                                                   # if not gready the unique right extension of reverse complement of x is x_

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
        self.SR.remove(x)
        if not is_palindromic(x):   self.SR.remove(get_reverse_sr(x))

        # 2
        if x != y:
            self.SR.remove(y)
            if not is_palindromic(y):   self.SR.remove(get_reverse_sr(y))

        # 3
        new=x+y[len_u:]
        t_new = tuple(new)
        if t_new not in self.abundances: 
            self.abundances[t_new]=(self.get_abundance(x)+self.get_abundance(y))          # TODO !!! how to deal with fusions for abundances
        else: 
            self.abundances[t_new]+=(self.get_abundance(x)+self.get_abundance(y))         # TODO !!! how to deal with fusions for abundances
        self.SR.sorted_add(new)
        if not is_palindromic(new): 
            rnew = get_reverse_sr(new)
            t_rnew = tuple(rnew)
            if t_rnew not in self.abundances: 
                self.abundances[t_rnew]=(self.get_abundance(x)+self.get_abundance(y))     # TODO !!! how to deal with fusions for abundances
            else: 
                self.abundances[t_rnew]+=(self.get_abundance(x)+self.get_abundance(y))    # TODO !!! how to deal with fusions for abundances
            
            
            self.SR.sorted_add(rnew)

        # we made a compaction, return 1.
        return 1
    
    def compaction(self,unitig_lengths,k,min_conflict_overlap):
        ''' Compaction of all sr in SR
        If a sr was not compacted, i will never be compacted.
        If it was compacted, maybe it will be re-compacted later on. However, no need to re-run the fusion on this read as
         - either I could have been compacted on the left and this was done before or this will be done latter or
         - it will be right extended later: the 'new' (the compacted sr) sr is higher in the lexicographic order than the original sr (as it is longer), thus its position is higher in the SR data structure, thus it will be seen again later.
        Note that this may be not true in parallel computations.
        '''

        checked=0
        compacted=0
        n = len(self.SR)
        for sr in self.SR.traverse():
            if checked%100==0:
                sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(self.SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\r")
            checked+=1
            witness = self.fusion(sr,unitig_lengths,k,min_conflict_overlap)
            if witness == 1: # a fusion was done
                compacted+=1
        sys.stderr.write("      Compacting, "+str(checked)+" checked. Size SR "+str(len(self.SR))+" %.2f"%(100*checked/n)+"%, "+str(compacted)+" couple of nodes compacted\n")

    def print_maximal_super_reads(self):
        '''print all maximal super reads as a flat format'''
        for sr in self.SR.traverse():
            if is_canonical(sr) or is_palindromic(sr):
                print(self.get_abundance(sr))
                if len(sr)==1:
                    print (str(sr[0])+";")
                else:
                    for unitig_id in sr:
                        print (str(unitig_id)+";", end="")
                    print ()
        
            
    def __str__(self):
        return "SR\n"+str(self.SR)+"\nabundances"+str(self.abundances)

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



def is_palindromic(x):
    ''' return true is a sr x is palindromic, eg [1,2,-2,-1]'''
    if len(x)%2 == 1: return False
    for i in range(0,int(len(x)/2)):
        if x[i] != -x[-i-1]: return False
    return True



def generate_SR(file_name):
    ''' Given an input file storing super reads, store them in the SR array'''
    ''' SR come with their abundance (one line over two) eg:
    62
    -201;-202;-640;-134;-358;-67;-68;-484;-48;
    SR are stored in sr and for each of them, it is associated to its abundance
    '''
    csr = counted_super_reads()
    with  open(file_name, 'r') as sr_file:
        while True:
            line = sr_file.readline()
            if not line: break
            abundance = int(line)
            
            line = sr_file.readline().rstrip()[:-1].split(';')
            sr=[]
            for unitig_id in line:
                sr_val=int(unitig_id)
                sr=sr+[sr_val]
            csr.add(sr, abundance)
    return csr
    
    
    





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



# def print_maximal_super_reads(SR):
#     '''print all maximal super reads as a flat format'''
#     for sr in SR.traverse():
#         if is_canonical(sr) or is_palindromic(sr):
#             if len(sr)==1:
#                 print (str(sr[0])+";")
#             else:
#                 for unitig_id in sr:
#                     print (str(unitig_id)+";", end="")
#                 print ()



def load_unitigs(file_name):
    unitigs=[]
    file=open(file_name,'r')
    for line in file:
        line=line.rstrip()
        if line[0]=='>': continue
        unitigs+=[line]
    return unitigs



def load_unitig_lengths(file_name):
    unitig_lengths=[]
    file=open(file_name,'r')
    for line in file:
        line=line.rstrip()
        if line[0]=='>': continue
        unitig_lengths+=[len(line)]
    return unitig_lengths



def get_len_ACGT_from_unitigs(msr,unitigs,size_overlap):
    lenACGT=0
    for unitig_id in msr:
        unitig_id = int(unitig_id)
        if unitig_id<0:                                         #the sequence is indicated as reverse complemented. Note that unitig ids start with 1, thus the -0 problem does not appear.
            unitig_id=-unitig_id
        lenACGT+=len(unitigs[unitig_id-1])
    lenACGT-=size_overlap*(len(msr) -1)                # remove the size of the overlaps
    return lenACGT



def get_len_ACGT(sr,unitig_lengths,size_overlap):
    ''' provides the cumulated length of the unitigs of a sr '''
    lenACGT = size_overlap
    for unitig_ids in sr:
        if unitig_ids>0:
            lenACGT+=unitig_lengths[unitig_ids-1]-size_overlap       # add the ACGT len of the corresponding unitig. -1 is due to the fact that unitigs start at zero while sr indices start at one.
        else:
            lenACGT+=unitig_lengths[-unitig_ids-1]-size_overlap

    return lenACGT






# def to_clean(SR,sr):
#     ''' A sr is a dead end if it has no successor or no predecessor '''
#     #~ if is_a_island(SR,sr):
#         #~ return True
#     if is_a_tip(SR,sr):
#         return True
#     return False


# def equivalent_context(SR,sr,succ,pred):
#     succ2=all_succ(SR,sr)
#     pred2=all_pred(SR,sr)
#     if(succ==succ2 and pred2==pred):
#     #~ if( all([z in succ2 for z in succ]) and all([z in pred2 for z in pred]) ):
#         return True


# def at_least_a_successor_with_equivalent_context(SR,sr,succ,pred,ref):
#     n=len(sr)
#     for len_u in range(n-1,0,-1):
#         Y=SR.get_lists_starting_with_given_prefix(sr[-len_u:])
#         for y in Y:
#             if(y!=ref):
#                 if equivalent_context(SR,y,succ,pred):
#                     return True
#     return False





# def clean_parallel_contigs(SR,sr):
#     succ=all_succ(SR,sr)
#     pred=all_pred(SR,sr)
#     if(len(pred)>0 and len(succ)>0):
#         for y in pred:
#             father=get_reverse_sr(y)
#             if(at_least_a_successor_with_equivalent_context(SR,father,succ,pred,sr)):
#                 #~ sys.stderr.write("FOUND IT OMG\n")
#                 SR.remove(sr)
#                 if not is_palindromic(sr): SR.remove(get_reverse_sr(sr))
#                 return





# def is_a_tip(SR,sr):
#     ''' A sr is a dead end if it has no successor or no predecessor '''
#     if is_a_dead_end(SR,sr):
#         if at_least_a_successor_with_multiple_predecessor(SR,sr):
#             return True
#         if at_least_a_successor_with_multiple_predecessor(SR,get_reverse_sr(sr)):
#             return True
#     return False



# def is_a_island(SR,sr):
#     ''' A sr is a dead end if it has no successor or no predecessor '''
#     if (not at_least_a_successor(SR,sr)) and  (not at_least_a_successor(SR,get_reverse_sr(sr))) :
#         return True
#     return False


# def is_a_dead_end(SR,sr):
#     ''' A sr is a dead end if it has no successor or no predecessor '''
#     if not at_least_a_successor(SR,sr):
#         return True # No predecesssor, this is a dead end
#     if not at_least_a_successor(SR,get_reverse_sr(sr)):
#         return True # No successor, this is a dead end
#     return False # Sucessor(s) and predecessor(s), this is not a dead end.
#
#


# def at_least_a_successor(SR,sr):
#     ''' Checks if sr as at least a successor. Return True in this case, else return False '''
#     n=len(sr)
#     for len_u in range(n-1,0,-1):
#         Y=SR.get_lists_starting_with_given_prefix(sr[-len_u:])
#         if (len(Y)>0):
#             return True
#     return False



# def at_least_a_successor_with_multiple_predecessor(SR,sr):
#     n=len(sr)
#     for len_u in range(n-1,0,-1):
#         Y=SR.get_lists_starting_with_given_prefix(sr[-len_u:])
#         for y in Y:
#             if multiple_successors(SR,get_reverse_sr(y)):
#                 return True
#     return False



# def all_succ(SR,sr):
#     res=[]
#     n=len(sr)
#     for len_u in range(n-1,0,-1):
#         Y=SR.get_lists_starting_with_given_prefix(sr[-len_u:])
#         for y in Y:
#             res.append(y)
#     res.sort()
#     return res



# def all_pred(SR,sr):
#     return all_succ(SR,get_reverse_sr(sr))


#Could be optimized by enumerating all 3 path check for a out if not continue etc
# def all_Qpaths(SR,sr,q):
#     if(q==0):
#         return all_succ(SR,sr)
#     res=[]
#     n=len(sr)
#     for len_u in range(n-1,0,-1):
#         Y=SR.get_lists_starting_with_given_prefix(sr[-len_u:])
#         for y in Y:
#             sons=all_Qpaths(SR,y,q-1)
#             for s in sons:
#                 s.insert(0,y)
#     res.sort()
#     return res


#possibly terrible performances
# def find_out_bulle(qpath,unitig_lengths,k):
#
#     if(len(qpath)==0):
#         return []
#     inter = [item for item in qpath[0] if (kc.get_len_ACGT(item,unitig_lengths,k) > 500)]
#     for X in qpath:
#         inter=list(set(inter),set(X))
#         if(inter.length()==0):
#             return []
#     return inter[0]



# def clean_complex_bulles(SR,sr,unitig_lengths,k,bulles_c):
#     qpath=all_Qpaths(SR,sr,bulles_c)
#     out=find_out_bulle(qpath,unitig_lengths,k)
#     if(not len(out)==0):
#         keep=qpath.pop()
#         for L in qpath:
#             for S in L:
#                 if not S in keep:
#                     SR.remove(get_reverse_sr(S))
#                     if not is_palindromic(sr): SR.remove(get_reverse_sr(sr))


# def multiple_successors(SR,sr):
#     n=len(sr)
#     for len_u in range(n-1,0,-1):
#         Y=SR.get_lists_starting_with_given_prefix(sr[-len_u:])
#         for y in Y:
#             if (at_least_a_successor_bis(SR,y)):
#                 return True
#     return False



# def at_least_a_successor_bis(SR,sr):
#     n=len(sr)
#     for len_u in range(n-1,0,-1):
#         Y=SR.get_lists_starting_with_given_prefix(sr[-len_u:])
#         for y in Y:
#             if at_least_a_successor(SR,y):
#                 return True
#     return False






def get_size_super_read_in_u(u,unitigs,k):
    '''compute the size of unitigs of a path stored in u (sum(size unitigs) - (|u|-1)*(k-1)'''
    cumulated_size_unitigs=0
    for unitig_id in u:
        if unitig_id<0: unitig_id=-unitig_id
        cumulated_size_unitigs+=len(unitigs[unitig_id-1])         # Ids starts at 1 (avoiding the -0=0 problem). Thus unitig i corresponds to unitigs[i-1]
    cumulated_size_unitigs-=(len(u)-1)*(k-1)                      # remove the size of the overlapping contigs
    return cumulated_size_unitigs



def get_msr_id(msr):
    ''' returns the id of a msr
    WARNING: here msr contains as last value its unique id.
    '''
    return int(msr[-1].split("_")[1])

def get_reverse_msr_id(msr,MSR):
    ''' returns the id of the reverse complement of msr
    1/ find the sequence of the msr_ in SR
    2/ grab its id
    WARNING: here msr contains as last value its unique id.
    '''
    #1/
    without_id_reverse_msr = get_reverse_sr(msr[:-1])                     # get the msr reverse complement
    Y=MSR.get_lists_starting_with_given_prefix(without_id_reverse_msr)        # find the reverse complement in the list.
    for y in Y:                                                               # should be of size >= 1. One of them is exactly 'without_id_reverse_msr' plus its id.
        if len(y) == len(without_id_reverse_msr)+1:                           # 'y' is 'without_id_reverse_msr' with its node id
            return get_msr_id(y)                                              # 2/
    return None                                                               # Should not happend

# SR=generate_SR("test.txt")
# SR.unique()
# SR=add_reverse_SR(SR)
# SR.unique()
# for sr in SR.traverse():
#     print (sr)
