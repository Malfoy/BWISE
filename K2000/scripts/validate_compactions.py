import sys


def get_rc(seq):
    rcseq=""
    for i in seq.split(";")[:-1]:
        if len(i)==0: continue
        
        if i[0]=='-':   i=i[1:]
        else:           i="-"+i
        rcseq=i+";"+rcseq
    return rcseq

compacted=sys.argv[1]
non_compacted=sys.argv[2]

compacted_text=""
for line in open(compacted): compacted_text+=line




for line in open(non_compacted):
    line=line.rstrip()
    if line not in compacted_text and get_rc(line) not in compacted_text:
        print (line+" KO")
    # else:
        # print (line+" ok")
