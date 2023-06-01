Origin of the TEs - repeat 2
================

# species identification

paper
<https://academic.oup.com/gbe/article/13/8/evab138/6307268#285072411>
Supplementary table Insect Genomes (Tables S1-S2).xlsx copy the table
into txt file: raw.txt

## Filter long reads and non Drosophila

``` bash
cat raw.txt |awk 'BEGIN{FS="\t"}$9=="long"'|grep -v "rosophila" > species-toinvestigateinfo.txt   
```

**Just retain the useful columns**

``` bash
cat species-toinvestigateinfo.txt|awk 'BEGIN{FS="\t";OFS="\t"}{print $7,$1,$2,$3,$5}' |perl -pe 's/ /./' > sampleids-toinvestigate.txt
#GCA_013387265.1    Coleoptera  Coccinellidae   Cryptolaemus    Cryptolaemus.montrouzieri
#GCA_011033045.1    Coleoptera  Coccinellidae   Harmonia    Harmonia.axyridis
#GCA_013421045.1    Coleoptera  Coccinellidae   Propylea    Propylea.japonica
#GCA_014170235.1    Coleoptera  Curculionidae   Listronotus Listronotus.bonariensis
#GCA_002938485.2    Coleoptera  Curculionidae   Sitophilus  Sitophilus.oryzae
#GCA_014611495.1    Coleoptera  Elateridae  Limonius    Limonius.californicus
#GCA_013368085.1    Coleoptera  Lampyridae  Abscondita  Abscondita.terminalis
#GCA_013368075.1    Coleoptera  Lampyridae  Lamprigera  Lamprigera.yunnana
#GCA_008802855.1    Coleoptera  Lampyridae  Photinus    Photinus.pyralis
#GCA_001937115.1    Coleoptera  Nitidulidae Aethina Aethina.tumida
#GCA_004143645.1    Coleoptera  Scarabaeidae    Protaetia   Protaetia.brevitarsis
#GCA_014905495.1    Coleoptera  Scarabaeidae    Trypoxylus  Trypoxylus.dichotomus
#GCA_001412225.1    Coleoptera  Silphidae   Nicrophorus Nicrophorus.vespilloides
#GCA_004115045.2    Collembola  Entomobryidae   Sinella Sinella.curviseta
#GCA_002217175.1    Collembola  Isotomidae  Folsomia    Folsomia.candida
#GCA_001718145.1    Collembola  Orchesellidae   Orchesella  Orchesella.cincta
#GCA_004302925.1    Diptera Calliphoridae   Cochliomyia Cochliomyia.hominivorax
#GCA_001735545.1    Diptera Calliphoridae   Phormia Phormia.regina
#GCA_000775305.1    Diptera Chironomidae    Belgica Belgica.antarctica
#GCA_001876365.2    Diptera Culicidae   Aedes   Aedes.albopictus
#GCA_002204515.1    Diptera Culicidae   Aedes   Aedes.aegypti
#GCA_013758885.1    Diptera Culicidae   Anopheles   Anopheles.albimanus
#GCA_001542645.1    Diptera Culicidae   Anopheles   Anopheles.gambiae
#GCA_003951495.1    Diptera Culicidae   Anopheles   Anopheles.funestus
#GCA_013141755.1    Diptera Culicidae   Anopheles   Anopheles.stephensi
#GCA_004136515.2    Diptera Culicidae   Anopheles   Anopheles.coluzzii
#GCA_002237135.1    Diptera Diopsidae   Teleopsis   Teleopsis.dalmanni
#GCA_003123925.1    Diptera Muscidae    Haematobia  Haematobia.irritans
#GCA_014843735.1    Diptera Muscidae    Musca   Musca.domestica
#GCA_014635995.1    Diptera Sarcophagidae   Sarcophaga  Sarcophaga.peregrina
#GCA_014529535.1    Diptera Sciaridae   Bradysia    Bradysia.coprophila
#GCA_001188975.4    Diptera Tephritidae Bactrocera  Bactrocera.oleae
#GCA_001854935.1    Hemiptera   Aleyrodidae Bemisia Bemisia.tabaci
#GCA_011764245.1    Hemiptera   Aleyrodidae Trialeurodes    Trialeurodes.vaporariorum
#GCA_014119065.1    Hemiptera   Anthocoridae    Orius   Orius.insidiosus
#GCA_009928515.1    Hemiptera   Aphididae   Aphis   Aphis.glycines
#GCA_003676215.3    Hemiptera   Aphididae   Rhopalosiphum   Rhopalosiphum.maidis
#GCA_008086715.1    Hemiptera   Aphididae   Sitobion    Sitobion.miscanthi
#GCA_003335185.2    Hemiptera   Delphacidae Laodelphax  Laodelphax.striatellus
#GCA_014356525.1    Hemiptera   Delphacidae Nilaparvata Nilaparvata.lugens
#GCA_000475195.1    Hemiptera   Liviidae    Diaphorina  Diaphorina.citri
#GCA_009739505.2    Hemiptera   Miridae Apolygus    Apolygus.lucorum
#GCA_003667255.1    Hemiptera   Pentatomidae    Euschistus  Euschistus.heros
#GCA_003261595.1    Hemiptera   Pseudococcidae  Maconellicoccus Maconellicoccus.hirsutus
#GCA_009761765.1    Hemiptera   Pseudococcidae  Phenacoccus Phenacoccus.solenopsis
#GCA_011037195.1    Hemiptera   Reduviidae  Triatoma    Triatoma.infestans
#GCA_003254395.2    Hymenoptera Apidae  Apis    Apis.mellifera
#GCA_009792835.1    Hymenoptera Apidae  Apis    Apis.dorsata
#GCA_003314205.1    Hymenoptera Apidae  Apis    Apis.mellifera mellifera
#GCA_011100585.1    Hymenoptera Apidae  Apis    Apis.cerana cerana
#GCA_013841205.1    Hymenoptera Apidae  Apis    Apis.mellifera caucasica
#GCA_013841245.1    Hymenoptera Apidae  Apis    Apis.mellifera carnica
#GCA_011952255.1    Hymenoptera Apidae  Bombus  Bombus.vosnesenskii
#GCA_011952205.1    Hymenoptera Apidae  Bombus  Bombus.bifarius
#GCA_011952275.1    Hymenoptera Apidae  Bombus  Bombus.vancouverensis nearcticus
#GCA_014905175.1    Hymenoptera Braconidae  Aphidius    Aphidius.gifuensis
#GCA_011426455.1    Hymenoptera Braconidae  Aphidius    Aphidius.ervi
#GCA_013357705.1    Hymenoptera Braconidae  Chelonus    Chelonus.insularis
#GCA_011426435.1    Hymenoptera Braconidae  Lysiphlebus Lysiphlebus.fabarum
#GCA_013123115.1    Hymenoptera Colletidae  Colletes    Colletes.gigas
#GCA_001855655.1    Hymenoptera Figitidae   Leptopilina Leptopilina.clavipes
#GCA_011634795.1    Hymenoptera Figitidae   Leptopilina Leptopilina.boulardi
#GCA_003227725.1    Hymenoptera Formicidae  Camponotus  Camponotus.floridanus
#GCA_009859135.1    Hymenoptera Formicidae  Formica Formica.selysi
#GCA_003227715.1    Hymenoptera Formicidae  Harpegnathos    Harpegnathos.saltator
#GCA_013373865.2    Hymenoptera Formicidae  Monomorium  Monomorium.pharaonis
#GCA_005281655.1    Hymenoptera Formicidae  Nylanderia  Nylanderia.fulva
#GCA_003672135.1    Hymenoptera Formicidae  Ooceraea    Ooceraea.biroi
#GCA_010367695.1    Hymenoptera Formicidae  Solenopsis  Solenopsis.invicta
#GCA_012274295.1    Hymenoptera Megachilidae    Osmia   Osmia.lignaria
#GCA_009193385.2    Hymenoptera Pteromalidae    Nasonia Nasonia.vitripennis
#GCA_012977825.2    Hymenoptera Pteromalidae    Pteromalus  Pteromalus.puparum
#GCA_010416925.1    Hymenoptera Vespidae    Polistes    Polistes.metricus
#GCA_010416935.1    Hymenoptera Vespidae    Polistes    Polistes.fuscatus
#GCA_014083535.1    Hymenoptera Vespidae    Vespa   Vespa.mandarinia
#GCA_014607495.2    Lepidoptera Carposinidae    Carposina   Carposina.sasakii
#GCA_004000445.1    Lepidoptera Crambidae   Chilo   Chilo.suppressalis
#GCA_014851415.1    Lepidoptera Crambidae   Cnaphalocrocis  Cnaphalocrocis.medinalis
#GCA_014595695.1    Lepidoptera Hesperiidae Epargyreus  Epargyreus.clarus clarus
#GCA_012273795.1    Lepidoptera Lasiocampidae   Dendrolimus Dendrolimus.punctatus
#GCA_002382865.1    Lepidoptera Noctuidae   Heliothis   Heliothis.virescens
#GCA_011316535.1    Lepidoptera Noctuidae   Spodoptera  Spodoptera.exigua
#GCA_011064685.1    Lepidoptera Noctuidae   Spodoptera  Spodoptera.frugiperda
#GCA_003590095.1    Lepidoptera Noctuidae   Trichoplusia    Trichoplusia.ni
#GCA_009667785.1    Lepidoptera Nymphalidae Maniola Maniola.jurtina
#GCA_009982905.1    Lepidoptera Pieridae    Colias  Colias.croceus
#GCA_005406025.1    Lepidoptera Psychidae   Eumeta  Eumeta.japonica
#GCA_004355975.1    Lepidoptera Pyralidae   Galleria    Galleria.mellonella
#GCA_014332785.1    Lepidoptera Saturniidae Antheraea   Antheraea.mylitta
#GCA_014132275.1    Lepidoptera Saturniidae Samia   Samia.ricini
#GCA_009982885.1    Lepidoptera Sphingidae  Hyles   Hyles.vespertilio
#GCA_014839805.1    Lepidoptera Sphingidae  Manduca Manduca.sexta
#GCA_003425675.2    Lepidoptera Tortricidae Cydia   Cydia.pomonella
#GCA_011170035.1    Orthoptera  Gryllidae   Teleogryllus    Teleogryllus.occipitalis
#GCA_003426905.1    Siphonaptera    Pulicidae   Ctenocephalides Ctenocephalides.felis
#GCA_012932325.1    Thysanoptera    Thripidae   Thrips  Thrips.palmi
#GCA_009617725.1    Trichoptera Hydropsychidae  Hydropsyche Hydropsyche.tenuis
#GCA_009617715.1    Trichoptera Polycentropodidae   Plectrocnemia   Plectrocnemia.conspersa
#GCA_008973525.1    Trichoptera Stenopsychidae  Stenopsyche Stenopsyche.tienmushanensis
```

# Download

Documentation
<https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/>

``` bash
cat sampleids-toinvestigate.txt|cut -f 1 >list
datasets download genome accession --inputfile list
unzip ncbi_dataset.zip 
mv **/*.fna fastafiles 
# muahaha -> nice one (abbreviate species name from "Bullshitus.giganticus" to B.giganticus)
cat sampleids-toinvestigate.txt|sort |awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$5}'|perl -pe 's/\s(.).+\./\t$1./'    
```

# Rename

``` bash
fasta-reader.py GCA_000475195.1_Diaci_psyllid_genome_assembly_version_1.1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/D.citri
fasta-reader.py GCA_000775305.1_ASM77530v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/B.antarctica
fasta-reader.py GCA_001188975.4_MU_Boleae_v2_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/B.oleae
fasta-reader.py GCA_001412225.1_Nicve_v1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/N.vespilloides
fasta-reader.py GCA_001542645.1_ASM154264v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.gambiae
fasta-reader.py GCA_001718145.1_ASM171814v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/O.cincta
fasta-reader.py GCA_001735545.1_ASM173554v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/P.regina
fasta-reader.py GCA_001854935.1_ASM185493v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/B.tabaci
fasta-reader.py GCA_001855655.1_ASM185565v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/L.clavipes
fasta-reader.py GCA_001876365.2_canu_80X_arrow2.2_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.albopictus
fasta-reader.py GCA_001937115.1_Atum_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.tumida
fasta-reader.py GCA_002204515.1_AaegL5.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.aegypti
fasta-reader.py GCA_002217175.1_ASM221717v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/F.candida
fasta-reader.py GCA_002237135.1_Tel_dalmanni_2A_v1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/T.dalmanni
fasta-reader.py GCA_002382865.1_K63_refined_pacbio_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/H.virescens
fasta-reader.py GCA_002938485.2_Soryzae_2.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/S.oryzae
fasta-reader.py GCA_003123925.1_Hi_v1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/H.irritans
fasta-reader.py GCA_003227715.1_Hsal_v8.5_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/H.saltator
fasta-reader.py GCA_003227725.1_Cflo_v7.5_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.floridanus
fasta-reader.py GCA_003254395.2_Amel_HAv3.1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.mellifera
fasta-reader.py GCA_003261595.1_MB_VBL1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/M.hirsutus
fasta-reader.py GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.mel.mel.
fasta-reader.py GCA_003335185.2_ASM333518v2_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/L.striatellus
fasta-reader.py GCA_003425675.2_Cpom.V2_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.pomonella
fasta-reader.py GCA_003426905.1_ASM342690v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.felis
fasta-reader.py GCA_003590095.1_tn1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/T.ni
fasta-reader.py GCA_003667255.1_E_heros_v1.0_sep_2018_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/E.heros
fasta-reader.py GCA_003672135.1_Obir_v5.4_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/O.biroi
fasta-reader.py GCA_003676215.3_ASM367621v3_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/R.maidis
fasta-reader.py GCA_003951495.1_AfunF3_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.funestus
fasta-reader.py GCA_004000445.1_ASM400044v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.suppressalis
fasta-reader.py GCA_004115045.2_ASM411504v2_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/S.curviseta
fasta-reader.py GCA_004136515.2_ASM413651v2_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.coluzzii
fasta-reader.py GCA_004143645.1_ASM414364v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/P.brevitarsis
fasta-reader.py GCA_004302925.1_ASM430292v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.hominivorax
fasta-reader.py GCA_004355975.1_Honeycombmoth_v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/G.mellonella
fasta-reader.py GCA_005281655.1_TAMU_Nfulva_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/N.fulva
fasta-reader.py GCA_005406025.1_Evar_2.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/E.japonica
fasta-reader.py GCA_008086715.1_ASM808671v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/S.miscanthi
fasta-reader.py GCA_008802855.1_Ppyr1.3_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/P.pyralis
fasta-reader.py GCA_008973525.1_ASM897352v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/S.tienmushanensis
fasta-reader.py GCA_009193385.2_Nvit_psr_1.1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/N.vitripennis
fasta-reader.py GCA_009617715.1_P_conspersa_v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/P.conspersa
fasta-reader.py GCA_009617725.1_HT1_wtdbg2_racon_pilon_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/H.tenuis
fasta-reader.py GCA_009667785.1_Mjurtina_v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/M.jurtina
fasta-reader.py GCA_009739505.2_ASM973950v2_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.lucorum
fasta-reader.py GCA_009761765.1_ASM976176v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/P.solenopsis
fasta-reader.py GCA_009792835.1_RUTG_Adors_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.dorsata
fasta-reader.py GCA_009859135.1_ASM985913v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/F.selysi
fasta-reader.py GCA_009928515.1_A_gly_BT4_v2.1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.glycines
fasta-reader.py GCA_009982885.1_ASM998288v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/H.vespertilio
fasta-reader.py GCA_009982905.1_SU_C_cro_Alba_ref_v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.croceus
fasta-reader.py GCA_010367695.1_ASM1036769v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/S.invicta
fasta-reader.py GCA_010416925.1_CU_Pmet_PB_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/P.metricus
fasta-reader.py GCA_010416935.1_CU_Pfus_HIC_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/P.fuscatus
fasta-reader.py GCA_011033045.1_ASM1103304v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/H.axyridis
fasta-reader.py GCA_011037195.1_UVM_Tinf_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/T.infestans
fasta-reader.py GCA_011064685.1_ZJU_Sfru_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/S.frugiperda
fasta-reader.py GCA_011100585.1_ASM1110058v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.cerana.cerana
fasta-reader.py GCA_011170035.1_Tocci_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/T.occipitalis
fasta-reader.py GCA_011316535.1_NJAU_Sexi_v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/S.exigua
fasta-reader.py GCA_011426435.1_ASM1142643v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/L.fabarum
fasta-reader.py GCA_011426455.1_ASM1142645v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.ervi
fasta-reader.py GCA_011634795.1_ASM1163479v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/L.boulardi
fasta-reader.py GCA_011764245.1_ASM1176424v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/T.vaporariorum
fasta-reader.py GCA_011952205.1_Bbif_JDL3187_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/B.bifarius
fasta-reader.py GCA_011952255.1_Bvos_JDL3184-5_v1.1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/B.vosnesenskii
fasta-reader.py GCA_011952275.1_Bvanc_JDL1245_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/B.vanc.nearcticus
fasta-reader.py GCA_012273795.1_ASM1227379v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/D.punctatus
fasta-reader.py GCA_012274295.1_USDA_OLig_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/O.lignaria
fasta-reader.py GCA_012932325.1_TpBJ-2018v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/T.palmi
fasta-reader.py GCA_012977825.2_ZJU_Ppup_chr_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/P.puparum
fasta-reader.py GCA_013123115.1_ASM1312311v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.gigas
fasta-reader.py GCA_013141755.1_UCI_ANSTEP_V1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.stephensi
fasta-reader.py GCA_013357705.1_ASM1335770v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.insularis
fasta-reader.py GCA_013368075.1_ASM1336807v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/L.yunnana
fasta-reader.py GCA_013368085.1_Ate_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.terminalis
fasta-reader.py GCA_013373865.2_ASM1337386v2_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/M.pharaonis
fasta-reader.py GCA_013387265.1_SYSU_Cmont_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.montrouzieri
fasta-reader.py GCA_013421045.1_ASM1342104v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/P.japonica
fasta-reader.py GCA_013758885.1_VT_AalbS3_pri_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.albimanus
fasta-reader.py GCA_013841205.1_ASM1384120v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.mel.caucasica
fasta-reader.py GCA_013841245.1_ASM1384124v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.mel.carnica
fasta-reader.py GCA_014083535.1_V.mandarinia_Nanaimo_p1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/V.mandarinia
fasta-reader.py GCA_014119065.1_ASM1411906v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/O.insidiosus
fasta-reader.py GCA_014132275.1_Sr_UT_JL_v1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/S.ricini
fasta-reader.py GCA_014170235.1_ASM1417023v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/L.bonariensis
fasta-reader.py GCA_014332785.1_AM_v1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.mylitta
fasta-reader.py GCA_014356525.1_ASM1435652v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/N.lugens
fasta-reader.py GCA_014529535.1_BU_Bcop_v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/B.coprophila
fasta-reader.py GCA_014595695.1_ASM1459569v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/E.clarus.clarus
fasta-reader.py GCA_014607495.2_ASM1460749v2_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.sasakii
fasta-reader.py GCA_014611495.1_UIdaho_Lcali_1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/L.californicus
fasta-reader.py GCA_014635995.1_ASM1463599v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/S.peregrina
fasta-reader.py GCA_014839805.1_JHU_Msex_v1.0_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/M.sexta
fasta-reader.py GCA_014843735.1_ASM1484373v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/M.domestica
fasta-reader.py GCA_014851415.1_ASM1485141v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/C.medinalis
fasta-reader.py GCA_014905175.1_ASM1490517v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/A.gifuensis
fasta-reader.py GCA_014905495.1_ASM1490549v1_genomic.fna | fasta-formatter.py --upper |fasta-writter.py > /Users/rokofler/analysis/dmel_TE_invasions/2023-05-otherspecies-ro/ncbi_dataset/data/rename/T.dichotomus
# executed with 
# cat torename|python renamer.py|zsh
```
