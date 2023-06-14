genome-proportion
================

# Assembly size

``` bash
for i in *.fasta; do echo $i; done                                                                                                                    
#Canton-S.fasta
#DGRP-732.fasta
#Iso1.fasta
#Oregon-R.fasta
#Pi2.fasta
#SZ244.fasta
#SZ45.fasta
for i in *.fasta; do samtools faidx $i; done 
for i in *.fasta.fai; do echo $i; cat $i|awk '{sum+=$2}END{print sum}'; done|paste - -                                                                  
#Canton-S.fasta.fai 149104628
#DGRP-732.fasta.fai 141550979
#Iso1.fasta.fai 143726002
#Oregon-R.fasta.fai 136257826
#Pi2.fasta.fai  167833895
#SZ244.fasta.fai    149637683
#SZ45.fasta.fai 174534457
```

# Proportion of TEs

Use defragmented annotation to avoid overlapping annotations and thus
counting a TE twice

``` bash
for i in *ori.out; do python ../info/rm-defragmenter.py --rm $i --fai ../info/teseqs.fasta.fai --dist 100 > ../raw-defrag/$i.def; done
```

**P-element**

``` bash
# note, $6 = divergence (we only want recent, i.e. non-diverged insertions)
for i in *ori.out.def; do echo $i; cat $i |awk '$1~/PPI251/'|awk '$6<2{s+=$7}END{print s}'; done |paste - -
Canton-S.fasta.ori.out.def  
DGRP-732.fasta.ori.out.def  40561
Iso1.fasta.ori.out.def  
Oregon-R.fasta.ori.out.def  78901
Pi2.fasta.ori.out.def   58027
SZ244.fasta.ori.out.def 
SZ45.fasta.ori.out.def  84921
```

**Hobo**

``` bash
# note, $6 = divergence (we only want recent, i.e. non-diverged insertions)
for i in *ori.out.def; do echo $i; cat $i |awk '$1~/DMHFL1/'|awk '$6<2{s+=$7}END{print s}'; done |paste - -
Canton-S.fasta.ori.out.def  
DGRP-732.fasta.ori.out.def  55253
Iso1.fasta.ori.out.def  32867
Oregon-R.fasta.ori.out.def  47898
Pi2.fasta.ori.out.def   188041
SZ244.fasta.ori.out.def 10069
SZ45.fasta.ori.out.def  17181
```

**I-element**

``` bash
# note, $6 = divergence (we only want recent, i.e. non-diverged insertions)
for i in *ori.out.def; do echo $i; cat $i |awk '$1~/DMIFACA/'|awk '$6<2{s+=$7}END{print s}'; done |paste - -
Canton-S.fasta.ori.out.def  34
DGRP-732.fasta.ori.out.def  122632
Iso1.fasta.ori.out.def  68158
Oregon-R.fasta.ori.out.def  103116
Pi2.fasta.ori.out.def   80656
SZ244.fasta.ori.out.def 102563
SZ45.fasta.ori.out.def  159571
```

**Tirant**

``` bash
# note, $6 = divergence (we only want recent, i.e. non-diverged insertions)
for i in *ori.out.def; do echo $i; cat $i |awk '$1~/TIRANT/'|awk '$6<2{s+=$7}END{print s}'; done |paste - -
Canton-S.fasta.ori.out.def  
DGRP-732.fasta.ori.out.def  84390
Iso1.fasta.ori.out.def  160922
Oregon-R.fasta.ori.out.def  5345
Pi2.fasta.ori.out.def   69154
SZ244.fasta.ori.out.def 20433
SZ45.fasta.ori.out.def  14532
```

**BLOOD**

``` bash
# note, $6 = divergence (we only want recent, i.e. non-diverged insertions)
for i in *ori.out.def; do echo $i; cat $i |awk '$1~/BLOOD/'|awk '$6<2{s+=$7}END{print s}'; done |paste - -
Canton-S.fasta.ori.out.def  260686
DGRP-732.fasta.ori.out.def  98543
Iso1.fasta.ori.out.def  218191
Oregon-R.fasta.ori.out.def  131216
Pi2.fasta.ori.out.def   259525
SZ244.fasta.ori.out.def 91783
SZ45.fasta.ori.out.def  217230
```

**OPUS**

``` bash
# note, $6 = divergence (we only want recent, i.e. non-diverged insertions)
for i in *ori.out.def; do echo $i; cat $i |awk '$1~/OPUS/'|awk '$6<2{s+=$7}END{print s}'; done |paste - -
Canton-S.fasta.ori.out.def  220601
DGRP-732.fasta.ori.out.def  167404
Iso1.fasta.ori.out.def  248259
Oregon-R.fasta.ori.out.def  278927
Pi2.fasta.ori.out.def   262659
SZ244.fasta.ori.out.def 21840
SZ45.fasta.ori.out.def  39518
```

**412**

``` bash
# note, $6 = divergence (we only want recent, i.e. non-diverged insertions)
for i in *ori.out.def; do echo $i; cat $i |awk '$1~/412/'|awk '$6<2{s+=$7}END{print s}'; done |paste - -
Canton-S.fasta.ori.out.def  249283
DGRP-732.fasta.ori.out.def  234147
Iso1.fasta.ori.out.def  241995
Oregon-R.fasta.ori.out.def  250280
Pi2.fasta.ori.out.def   312390
SZ244.fasta.ori.out.def 69769
SZ45.fasta.ori.out.def  168469
```

## tmp - old

``` bash
cat Pi2.fasta.ori.out.def|awk '$1~/PPI251|TIRANT|OPUS|412|BLOOD|DMIFACA|DMHFL1/'|awk '$6<2{s+=$7}END{print s}'
```
