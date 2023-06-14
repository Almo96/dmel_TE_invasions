#!/usr/bin/env python
import os
import sys
import re
import argparse
import random
import collections




parser = argparse.ArgumentParser(description="""           
Description
-----------
    This script simulates single-end reads from the population genome""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 3+

Authors
-------
    Robert Kofler 
""")



parser.add_argument("--rm", type=str, required=True, dest="rm", default=None, help="the repeatmasker file")
args = parser.parse_args()

teh=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
temax=collections.defaultdict(lambda:0)
speciesset=set([])



for line in open(args.rm):
	line=line.rstrip().lstrip()
	a=re.split("\s+",line)
	#print(a)
	#  0		1		2		3		4			5		6			7			8	9			10				11		12		13		14		
	#['277', '25.72', '0.00', '1.43', 'contig_1', '25830', '25900', '(34454495)', 'C', 'BAGGINS', 'Unspecified', '(1149)', '4304', '4235', 'C.costata.fasta.ori.out']
	#['256', '29.38', '2.03', '3.42', 'contig_1', '32698', '32845', '(34447550)', '+', 'INE1', 'Unspecified', '29', '174', '(437)', 'C.costata.fasta.ori.out']
	score,te,species=int(a[0]),a[9],a[14]
	speciesset.add(species)
	if score> temax[te]:
		temax[te]=score
	if score > teh[te][species]:
		teh[te][species] =score

for te in temax.keys():
	dmelscore=temax[te]
	if dmelscore<1:
		continue
	for species in speciesset:
		relscore=float(teh[te][species])/float(dmelscore)
		print("{0}\t{1}\t{2}".format(te,species,relscore))
