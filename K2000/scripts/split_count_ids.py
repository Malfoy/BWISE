import sys

file=open(sys.argv[0])
for line in file:
	for id in line.rstrip().split(";")[:-1]:
		if id[0]=='-':
			print(id[1:])
		else:
			print (id)
