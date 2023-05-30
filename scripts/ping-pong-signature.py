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


def get_end(start,cig):
     
     result = re.findall(prog,cig)
     alignmentleng=0
     for count,cigchar in result:
          count=int(count)
          if cigchar=="M" or cigchar=="D":
               alignmentleng+=count
     
     end=start+alignmentleng-1
     return end

def load_5p_position(file,minmq,maxmm,gtag,minlen,maxlen):
     ps= collections.defaultdict(lambda:collections.defaultdict(lambda:0))
     pas=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
     tes=set([])
     for line in file:
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
          tmp=a[11]
          b=tmp.split(" ")
          for bt in b:
               if bt.startswith("NM:i:"):
                    mm=int(bt[5:])
          if(mm>maxmm):
               continue  
          
          
          ref=a[2]
          readlen=len(a[9])


          if ref.endswith(gtag):
               le=len(gtag)
               teseq=ref[:-le]
               
               if readlen<minlen or readlen>maxlen:
                    continue
               tes.add(teseq)
               start=int(a[3])
               if flag& 0x10:
                    # reverse complement if flag 0x10 is set
                    end=get_end(start,a[5]) # for reverse complements get the end
                    pas[teseq][end]+=1
               else:
                    ps[teseq][start]+=1
                         

     return tes,ps,pas


# already for a given TE
def get_s_sig(sense,anti,siglen,minol):
     sig=[0.0]*siglen
     
     if len(sense)<1 or len(anti)<1:
          return None
     
     ssum=0.0

     for pos,count in sense.items():
          ssum+=float(count)
          
     if ssum<1:
          return None
     olc=0.0
     
     for i,rawcount in sense.items():
         
          weight=0.0
          weight=float(rawcount)/ssum
               
          #123456789
          #|
               
          for k in range(0,siglen):
               pos=i+k
               if pos in anti:
                    toadd=anti[pos]*weight
                    sig[k]+=toadd
                    olc+=max([anti[pos],rawcount])
     
     if olc<minol:
          return None

     return normalizecounts(sig)




def get_as_sig(sense,anti,siglen,minol):
     sig=[0.0]*siglen
     
     if len(sense)<1 or len(anti)<1:
          return None
     
     assum=0.0
     for pos,count in anti.items():
          assum+=float(count)    
     if assum<1:
          return None
     olc=0.0
     
     for i,rawcount in anti.items():
         
          weight=0.0
          weight=float(rawcount)/assum
               
          #     * len3
          #01234567
          #   ^  3=5-3+1
          #     
               
          sensestart=i-siglen+1  
          for k in range(0,siglen):
               pos=sensestart+k
               if pos in sense:
                    toadd=sense[pos]*weight
                    sig[k]+=toadd
                    olc+=max([sense[pos],rawcount])
     
     if olc<minol:
          return None

     return normalizecounts(sig)


def normalizecounts(sig):
     sigsum=float(sum(sig))
     for i in range(0,len(sig)):
          if sigsum>0.0:
               sig[i]/=sigsum
     return sig

def load_exlcudesites(ip):
     toret=[]
     if ip=="":
          pass
     elif("," in ip):
          a=ip.split(",")
          toret=[int(i) for i in a]
     else:
          toret.append(int(ip))
     return set(toret)
          
     

parser = argparse.ArgumentParser(description="""           
Description
-----------
Summary statistics
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
miRNA: 21-23nt
piRNA: 23-28nt


Authors
-------
    Robert Kofler
""")
parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="A sam file")
parser.add_argument("--min-mq", type=int, required=False, dest="minmq", default=1, help="min mapping quality")
parser.add_argument("--max-mm", type=int, required=False, dest="maxmm", default=2, help="max mismatches")
parser.add_argument("--min-len", type=int, required=False, dest="minlen", default=23, help="min length")
parser.add_argument("--max-len", type=int, required=False, dest="maxlen", default=29, help="max length")
parser.add_argument("--sig-len", type=int, required=False, dest="siglen", default=19, help="signature length")
parser.add_argument("--min-ol", type=int, required=False, dest="minol", default=5, help="minimum number of overlapping ping-pong pairs")
parser.add_argument("--sample-id", type=str, required=True, dest="sid", default=10, help="the sample id")
parser.add_argument("--group-tag", type=str, required=False, dest="group", default="_te", help="the sample group, e.g. _te, _mRNA, _tRNA")
parser.add_argument("--exclude-sites", type=str, required=False, dest="exclude", default="", help="a comma separated list of sites to exclude from the pelement")

args = parser.parse_args()
minmq=args.minmq
maxmm=args.maxmm
exsite=load_exlcudesites(args.exclude)
minlen,maxlen = args.minlen, args.maxlen
tes,ps,pas = load_5p_position(args.sam,args.minmq,args.maxmm,args.group,minlen,maxlen)


# Mask the P-element sites
# TODO: generalized for all other TE families
for i in exsite:
     ps["PPI251"][i]=0
     pas["PPI251"][i]=0

sid=args.sid
for te in tes:
     sense=ps[te]
     anti=pas[te]
     
     #if te =="PPI251":
     #     print sense
     #     print anti
     asar=get_as_sig(sense,anti,args.siglen,args.minol)
     sar=get_s_sig(sense,anti,args.siglen,args.minol)
     
     if asar is not None:
          for i in range(0,len(asar)):
               print("{0}\t{1}\t{2}\t{3}\t{4}".format(sid,te,"as",i+1,asar[i]))

     if sar is not None:
          for i in range(0,len(sar)):
               print("{0}\t{1}\t{2}\t{3}\t{4}".format(sid,te,"s",i+1,sar[i]))
