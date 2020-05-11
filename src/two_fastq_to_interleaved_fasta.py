import sys
import gzip
"""
Given two fq files of equal number of reads, prints a unique fasta file interleaving reads from input fastq files
"""

__author__      = "Pierre Peterlongo"
__email__       = "pierre.peterlongo@inria.fr"


def comp(letter):
    if letter=='A':
        return 'T'
    if letter=='T' :
        return 'A'
    if letter=='G' :
        return 'C'
    if letter=='C' :
        return 'G'
    if letter=='a':
        return 't'
    if letter=='t' :
        return 'a'
    if letter=='g' :
        return 'c'
    if letter=='c' :
        return 'g'
    return letter
    

def revcomp(seq):
    res=""
    for letter in range(len(seq)-1,-1,-1):
        res=res+str(comp(seq[letter]))
    return res

def deinterleave (fileNameA: str, fileNameB: str):
    # fileA = open(fileNameA,"r")
    try:
        fileA = gzip.open(fileNameA, 'rt')
        fileA.read(1)
    except OSError:
        fileA = open(fileNameA,"r")
    if not fileA:
        raise IOError(f"File {fileNameA} cannot be opened")
    fileA.seek(0)

    
    try:
        fileB = gzip.open(fileNameB, 'rt')
        fileB.read(1)
    except OSError:
        fileB = open(fileNameB,"r")
    if not fileB:
        raise IOError(f"File {fileNameB} cannot be opened")
    fileB.seek(0)

    while True: 
        lineA=fileA.readline()
        lineB=fileB.readline()
        if not lineA and lineB: 
                raise ValueError(f"Files {fileNameA} and {fileNameB} must contain the same number of lines")
        if lineA and not lineB: 
                raise ValueError(f"Files {fileNameA} and {fileNameB} must contain the same number of lines")
        if not lineA: break
        lineA = lineA.strip()
        lineB = lineB.strip()
        sequenceA = fileA.readline().strip()
        sequenceB = fileB.readline().strip()
        fileA.readline().strip()      # don't care about the second comment line
        fileA.readline().strip()      # don't care about the quality
        fileB.readline().strip()      # don't care about the second comment line
        fileB.readline().strip()      # don't care about the quality
        print(f">{lineA[1:]}\n{sequenceA}\n>{lineB[1:]}\n{revcomp(sequenceB)}")
    fileA.close()
    fileB.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:

        print(f"Usage: python {sys.argv[0]} file1.fq[.gz] file2.fq[.gz]\n\t Interleave the two fastq files (possibly gzipped).\n\t Reads from file2 are reverse complemented \n\t Prints the result on the standard output. ")
        exit()
    deinterleave(sys.argv[1], sys.argv[2])