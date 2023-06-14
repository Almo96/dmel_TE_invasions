#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import math
import collections




class TEInsertion:
	tecounter=1
	def __init__(self,fam,chrm,start,end,strand,div,qstart,qend,sid):
		self.fam=fam
		self.chrm=chrm
		self.start=start
		self.strand=strand
		self.end=end
		self.div=div
		self.qstart=qstart
		self.qend=qend
		self.sid=sid
		self.counter=TEInsertion.tecounter
		TEInsertion.tecounter+=1
	
	def __str__(self):
		# chr1 20 25 forward 1 +
		# chr1 20 25 reverse 1 -
		name=self.fam+"_"+self.sid+"_"+str(self.counter)
		#lengsq=str(self.qend-self.qstart)
		lengsref=str(self.end-self.start)

		topr=[self.chrm, str(self.start), str(self.end),name ,lengsref,  self.strand]
		return "\t".join(topr)


"""
class TEBuilder:
	def __init__(self,length,fam,strand,chrm,start):
		self._seq=[0]*(length+1)
		self._con=[]
		self.maxend=0
		self._len=length
		self._fam=fam
		self._strand=strand
		self._chrm=chrm
		self._start=start

	def add(self,te):
		assert(te.fam==self._fam)
		assert(te.strand==self._strand)
		assert(te.chrm==self._chrm)
		start,end,div=te.qstart,te.qend,te.div
		self._con.append((start,end,div))
		for i in range(start,end+1):
			self._seq[i]=1

		if te.end > self.maxend:
			self.maxend=te.end
	
	def length(self):
		return sum(self._seq)

	def avdiv(self):
		score=0.0
		tol=0
		for s,e,d in self._con:
			ems=e-s+1
			tol+=ems
			score+=ems*d
		av=score/float(tol)
		return av
	
	def __str__(self):
		topr=[self._fam, self._chrm, self._strand, str(self._start), str(self.maxend),str(self.avdiv()),str(self.length()),str(self._len)]
		return "\t".join(topr)
"""


def read_rm(rmfile,sid):
	reli=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[])))
	for l in open(rmfile):
		l=l.strip()
		a=re.split("\s+",l)
		# 0			1		2		3		4					5	6		7		8		9			10			11			12		13
		# ['997', '4.88', '0.00', '0.00', '211000022278031', '345', '467', '(554)', 'C', 'DM06920', 'Unspecified', '(4358)', '1725', '1603']
		score,div,chrm,start,end,strand,fam = int(a[0]),float(a[1]),a[4],int(a[5]),int(a[6]),a[8],a[9]
		raw1,raw2,raw3=a[11],a[12],a[13]
		qstart,qend=0,0
		if(strand=="C"):
			qstart,qend=int(raw3),int(raw2)
			strand="-"
		elif(strand=="+"):
			qstart,qend=int(raw1),int(raw2)
		else:
			raise Exception("Invalid")
		#  26532 0.00 0.31 0.31             2L_RaGOO   3975229   3978135 (22168013) +               PPI251     Unspecified       1    2907     (0)    
 		#   1117 0.82 0.00 0.00             2L_RaGOO   4095206   4095327 (22050821) +               PPI251     Unspecified       1     122  (2785)  
		#  26711 0.03 0.24 0.17             2L_RaGOO  25681458  25684362 (461786) C               PPI251     Unspecified     (0)    2907       1    
		#   4617 0.38 0.57 0.00             2L_RaGOO  24790272  24790799 (1355349) 	C               PPI251     Unspecified     (0)    2907    2377    
  		#	3610 0.00 0.51 0.00             2L_RaGOO  24790789  24791182 (1354966) 	C               PPI251     Unspecified  (2493)     414      19    
		assert(end>start)
		assert(qend>qstart)
		te=TEInsertion(fam,chrm,start,end,strand,div,qstart,qend,sid)
		reli[fam][chrm][strand].append(te)
	return reli

 

parser = argparse.ArgumentParser(description="""           
Description
-----------
Summary statistics
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""

Authors
-------
    Robert Kofler
""")
parser.add_argument("--rm", type=str, required=True, dest="rm", default=None, help="repeat masker file")
parser.add_argument("--sid", type=str, required=True, dest="sid", default=None, help="sample ID")
parser.add_argument("--dist", type=int, required=False, dest="dist", default=10, help="maximum distance between consecutive matches of the same TE")
#parser.add_argument("--dist", type=int, required=True, dest="dist", default=None, help="distance")
args = parser.parse_args()


rm=read_rm(args.rm,args.sid)
maxdist=args.dist


for fam, tmp2 in rm.items():
	for chrm, tmp1 in tmp2.items():
		for strand, telist in tmp1.items():
			al=sorted(telist,key=lambda i:i.start)
			active=al.pop(0)
			for te in al:
				if active.end+maxdist >= te.start:
					active.end=te.end
					if te.qend > active.qend:
						active.qend=te.qend
					if te.qstart < active.qstart:
						active.qstart=te.qstart
				else:
					print(str(active))
					active=te
			print(str(active))






