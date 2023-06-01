#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import math
import collections
import fileinput
prog = re.compile(r"(\d+)([MISDHN])")




parser = argparse.ArgumentParser(description="""           
Description
-----------
Summary statistics
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Note: normalization in this script is to a million piRNAs (normally it would be a million miRNAs)
piRNAs are defined as sequences
- mapping to a TE (reference name ending with *_te)
- and between the --pi-min and --pi-max

miRNA: 21-23nt
piRNA: 23-28nt


Authors
-------
    Robert Kofler
""")
parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="A sam file")
parser.add_argument("--min-mq", type=int, required=False, dest="minmq", default=0, help="min mapping quality")
parser.add_argument("--max-mm", type=int, required=False, dest="maxmm", default=2, help="max mismatches")
parser.add_argument("--pi-min", type=int, required=False, dest="pistart", default=23, help="minimum size of piRNAs")
parser.add_argument("--pi-max", type=int, required=False, dest="piend", default=29, help="maximum size of piRNAs")
parser.add_argument("--sample-id", type=str, required=True, dest="sid", default=10, help="the sample id")
args = parser.parse_args()
minmq=args.minmq
maxmm=args.maxmm


mistart=21
miend=23
pistart=args.pistart
piend=args.piend

tecount=collections.defaultdict(lambda:0)
tesum=0
mirnacount=0
trnacount=0
rrnacount=0
pirnacount=0
pirnanorm=0


ps=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
pas=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
tes=set([])

reo=re.compile(r"NM:i:(\d+)")

def get_end(start,cig):
     
     result = re.findall(prog,cig)
     alignmentleng=0
     for count,cigchar in result:
          count=int(count)
          if cigchar=="M" or cigchar=="D":
               alignmentleng+=count
     
     end=start+alignmentleng-1
     return end
          


for line in args.sam:
     """
0         1         2              3    4         5    6         7      8            9                        10                  11
r1	16	M14653_te	172	70	23M	*	0	0	ATGTCGAGTTTCGTGCCGAATAA	FFFFFFFFFFFFFFFFFFBBBBB	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:23
r2	0	M14653_te	240	70	27M	*	0	0	AACAGCTGCGGAATCGCACCGAATGCT	BBBBBFFFFFBFFFFFFFFFFFFFFFF	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:27
     """
     a=line.rstrip("\n").split("\t")
     
     # discard unmapped
     flag=int(a[1])
     if flag & 0x004 > 0:
          continue
     
     # discard low mapping quality
     mq=int(a[4])
     if mq< minmq:
          continue
     
     # discard mismatch
     mm=0
     mo=re.search(reo,line)
     if mo is None:
          continue
     
     mm=int(mo.group(1))
     if(mm>maxmm):
          continue  
     ref=a[2]
     readlen=len(a[9])
     if ref.endswith("_te"):
          if readlen>=23 and readlen<=29:
               pirnanorm+=1  # normalization always to the ones between 23 and 29
               # therefore no continue here and hard coding of sm
               
          if readlen<pistart or readlen>piend:
               continue # ignore everything smaller than pistart and larger than piend
          teseq=ref[:-3]
          tes.add(teseq)
          
          start=int(a[3])
          if flag& 0x10:
               # reverse complement if flag 0x10 is set
               end=get_end(start,a[5]) # for reverse complements get the end
               pas[teseq][end]+=1
          else:
               ps[teseq][start]+=1
                    
     elif ref.endswith("_miRNA"):
          if readlen<mistart or readlen>miend:
               continue
          mirnacount+=1
     elif ref.endswith("_rRNA") or  ref.endswith("_rRNA;"):
          rrnacount+=1
     elif ref.endswith("_tRNA"):
          trnacount+=1
     elif  ref.endswith("_snoRNA;") or ref.endswith("_snoRNA") or ref.endswith("_snRNA;") or ref.endswith("_snRNA") or ref.endswith("_mRNA"):
          pass
     else:
          raise Exception("Unknown sequence end "+ ref)




sid=args.sid
for te in tes:
     for start,count in ps[te].items():
          normcount=float(count)/float(pirnanorm)
          normcount*=1000000.0
          print("{0}\t{1}\t{2}\t{3}".format(sid,te,start,normcount))
          
for te in tes:
     for start,count in pas[te].items():
          normcount=float(count)/float(pirnanorm)
          normcount*=1000000.0
          print("{0}\t{1}\t{2}\t{3}".format(sid,te,start,-normcount))





