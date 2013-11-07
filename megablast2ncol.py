#!/bin/usr/python

#transform megablast output in ncol format and filter out hits with a match accuracy and a matchlength below the given threshold
#The script reads the output of megablast (-D 3 option for output format).
#It takes a lot of memory if the megablast output file is big.

import sys
import string
from operator import itemgetter

def megablast2ncol(fastaFile, alignFile, outFile, matchaccuracy, matchlength):

	fastaF = open(fastaFile, "r")
	alignF = open(alignFile, "r")
	outF = open(outFile, "w")

	dicFasta = {}
	seqA = set()
	for line in fastaF:
		if(line[0] == ">"):
			seqname = string.split(string.split(line, ">")[1],"\n")[0]
			seqname = seqname.split()[0]
			seqA.add(seqname)
		else:
			seq = string.split(line,"\n")[0]
			try:
				dicFasta[seqname]
			except:
				dicFasta[seqname] = len(seq)
	fastaF.close()

	#read output of megablast (-D 3 option for output format has been used)
	dicEdges = {}
	edgesList = []

	for line in alignF:
		if (line[0] != "#"):
			l = line.split()
			q = l[0]
			t = l[1]
			if (t != q):
				iden = float(l[2])/100
				align_len = int(l[3])
				mismatches = int(l[4])
				bitScore = int(float(l[11])) #round off to int
				q_len = dicFasta[q]
				t_len = dicFasta[t]
				
				coverage = align_len/float(q_len)
				
				if (iden >= matchaccuracy) and (coverage >= matchlength):
					try:
						dicEdges[q,t]
						dicEdges[t,q]
					except:
						dicEdges[q,t] = bitScore
						dicEdges[t,q] = bitScore
						edgesList.append([q,t,bitScore])
	alignF.close()
	del(dicFasta)
	del(dicEdges)

	edgesList_sorted = sorted(edgesList, key=itemgetter(2), reverse=True)
	del(edgesList)
	
	totNode = set()
	for edge in edgesList_sorted:
		outF.write(edge[0] + "\t" + edge[1] + "\t" + str(edge[2]) + "\n")
		totNode.add(edge[0])
		totNode.add(edge[1])
	outF.close()

	del(edgesList_sorted)

	print len(totNode)

fastaFile = sys.argv[1]
alignFile = sys.argv[2]
outFile = sys.argv[3]
matchaccuracy = float(sys.argv[4])
matchlength = float(sys.argv[5])

megablast2ncol(fastaFile, alignFile, outFile, matchaccuracy, matchlength)


