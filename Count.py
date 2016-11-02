#!/usr/bin/env python
import subprocess,sys,struct,os,time,shutil,random

sequence = ""
nuc=0
seq=0
with open(sys.argv[1]) as infile:  # do not load the whole file, one line at time
	for line in infile:
		string = line.rstrip()
		if ">" not in string:
			nuc+=len(line)
			seq+=1
			#~ print len(line)
		#~ else:
			#~ print string
	print "Nucleotide number"
	print nuc
	print "Sequence number"
	print seq
	print "Mean length of sequences"
	print nuc/seq
