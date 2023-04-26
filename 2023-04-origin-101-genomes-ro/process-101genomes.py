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

specieslist = [ "D.melanogaster.fasta.ori.out","D.simulans.fasta.ori.out","D.mauritiana.fasta.ori.out","D.sechellia.fasta.ori.out",
"D.yakuba.fasta.ori.out","D.teissieri.2733.fasta.ori.out","D.teissieri.CT02.fasta.ori.out","D.erecta.fasta.ori.out","D.ananassae.fasta.ori.out",
"C.costata.fasta.ori.out","D.acanthoptera.fasta.ori.out","D.ambigua.fasta.ori.out","D.americana.fasta.ori.out","D.arawakana.fasta.ori.out",
"D.biarmipes.fasta.ori.out","D.bipectinata.fasta.ori.out","D.cardini.fasta.ori.out","D.carrolli.fasta.ori.out","D.dunni.fasta.ori.out","D.elegans.fasta.ori.out","D.equinoxialis.fasta.ori.out",
"D.ercepeae.fasta.ori.out","D.eugracilis.fasta.ori.out","D.ezoana.masurca.fasta.ori.out","D.ficusphila.fasta.ori.out","D.funebris.fasta.ori.out","D.fuyamai.fasta.ori.out",
"D.grimshawi.fasta.ori.out","D.immigrans.fk05.fasta.ori.out","D.immigrans.kari17.fasta.ori.out","D.insularis.fasta.ori.out","D.jambulina.fasta.ori.out","D.kikkawai.fasta.ori.out",
"D.kurseongensis.fasta.ori.out","D.littoralis.fasta.ori.out","D.m.malerkotliana.fasta.ori.out","D.m.pallens.fasta.ori.out",
"D.mojavensis.fasta.ori.out","D.murphyi.fasta.ori.out","D.nebulosa.fasta.ori.out","D.neocordata.fasta.ori.out","D.obscura.fasta.ori.out","D.oshimai.fasta.ori.out","D.p.nigrens.fasta.ori.out",
"D.p.pseudoananassae.fasta.ori.out","D.parabipectinata.fasta.ori.out","D.paulistorum.L06.fasta.ori.out","D.paulistorum.L12.fasta.ori.out","D.persimilis.fasta.ori.out","D.prosaltans.fasta.ori.out",
"D.pruinosa.fasta.ori.out","D.pseudoobscura.fasta.ori.out","D.quadrilineata.fasta.ori.out","D.repleta.fasta.ori.out","D.repletoides.fasta.ori.out","D.rhopaloa.fasta.ori.out","D.rufa.fasta.ori.out",
"D.saltans.fasta.ori.out","D.sp.mushsaotome.fasta.ori.out","D.spaffchauvacae.fasta.ori.out","D.sproati.fasta.ori.out","D.sturtevanti.fasta.ori.out",
"D.subobscura.fasta.ori.out","D.subpulchrella.L1.fasta.ori.out","D.sucinea.fasta.ori.out","D.takahashii.fasta.ori.out",
"D.triauraria.fasta.ori.out","D.tristis.fasta.ori.out","D.tropicalis.fasta.ori.out","D.varians.fasta.ori.out","D.virilis.fasta.ori.out","D.wassermani.fasta.ori.out","D.willistoni.00.fasta.ori.out",
"D.willistoni.17.fasta.ori.out","Drosophila.willistoni.17.fasta.ori.out","L.clarofinis.fasta.ori.out","L.collinella.fasta.ori.out","L.magnipectinata.fasta.ori.out",
"L.mommai.fasta.ori.out","L.stackelbergi.fasta.ori.out","L.varia.fasta.ori.out","S.graminum.fasta.ori.out","S.hsui.fasta.ori.out","S.montana.fasta.ori.out","S.pallida.fasta.ori.out",
"Z.africanus.fasta.ori.out","Z.camerounensis.fasta.ori.out","Z.capensis.fasta.ori.out","Z.davidi.fasta.ori.out","Z.gabonicus.fasta.ori.out","Z.ghesquierei.fasta.ori.out",
"Z.indianus.16GNV01.fasta.ori.out","Z.indianus.BS02.fasta.ori.out","Z.indianus.CDD18.fasta.ori.out","Z.indianus.RCR04.fasta.ori.out","Z.inermis.fasta.ori.out","Z.kolodkinae.fasta.ori.out",
"Z.lachaisei.fasta.ori.out","Z.nigranus.fasta.ori.out","Z.ornatus.fasta.ori.out","Z.taronus.fasta.ori.out","Z.tsacasi.car7-4.fasta.ori.out","Z.tsacasi.fasta.ori.out","Z.vittiger.fasta.ori.out"]


for line in open(args.rm):
	line=line.rstrip().lstrip()
	a=re.split("\s+",line)
	#  0		1		2		3		4			5		6			7			8	9			10				11		12		13		14		
	#['277', '25.72', '0.00', '1.43', 'contig_1', '25830', '25900', '(34454495)', 'C', 'BAGGINS', 'Unspecified', '(1149)', '4304', '4235', 'C.costata.fasta.ori.out']
	#['256', '29.38', '2.03', '3.42', 'contig_1', '32698', '32845', '(34447550)', '+', 'INE1', 'Unspecified', '29', '174', '(437)', 'C.costata.fasta.ori.out']
	score,te,species=int(a[0]),a[9],a[14]
	if score> temax[te]:
		temax[te]=score
	if score > teh[te][species]:
		teh[te][species] =score

for te in temax.keys():
	dmelscore=temax[te]
	if dmelscore<1:
		continue
	for species in specieslist:
		relscore=float(teh[te][species])/float(dmelscore)
		print("{0}\t{1}\t{2}".format(te,species,relscore))
