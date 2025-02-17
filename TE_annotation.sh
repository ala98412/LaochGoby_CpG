#De novo
/home/why/Tools/RepeatMasker_Tools/RepeatModeler-2.0.2a/BuildDatabase -name RA_polished_DB -engine rmblast RA.final_polished.fasta
/home/why/Tools/RepeatMasker_Tools/RepeatModeler-2.0.2a/RepeatModeler -engine rmblast -database RA_polished_DB -pa 10 &> De_novo_RA_polished.run.out
/home/why/Tools/RepeatMasker_Tools/RepeatMasker/RepeatMasker -s -pa 10 -gff -engine rmblast -lib RM_587686.ThuAug111114132022/consensi.fa.classified -xsmall RA.final_polished.fasta

#Fish TEDB (PM)
#Periophthalmus_magnuspinnatus_transposable_elements (PM_TEDB)
/home/why/Tools/RepeatMasker_Tools/RepeatMasker/RepeatMasker -s -pa 10 -gff -engine rmblast -lib Periophthalmus_magnuspinnatus_transposable_elements -xsmall RA.final_polished.fasta

#RepBase + Dfam
#/home/why/Tools/RepeatMasker_Tools/RepeatMasker/famdb.py -i /home/why/Tools/RepeatMasker_Tools/RepeatMasker/Libraries/RepeatMaskerLib.h5 families --format fasta_name --ancestors --descendants 'Percomorphaceae' --include-class-in-name > Percomorphaceae_TE_lib.fa
/home/why/Tools/RepeatMasker_Tools/RepeatMasker/RepeatMasker -s -pa 10 -gff -engine rmblast -lib Percomorphaceae_TE_lib.fa -xsmall RA.final_polished.fasta


#LTR
#!/bin/bash

RepeatModeler=/home/why/Tools/RepeatMasker_Tools/RepeatModeler-2.0.2a/RepeatModeler
BuildDatabase=/home/why/Tools/RepeatMasker_Tools/RepeatModeler-2.0.2a/BuildDatabase
RepeatClassifier=/home/why/Tools/RepeatMasker_Tools/RepeatModeler-2.0.2a/RepeatClassifier
repeatmasker=/home/why/Tools/RepeatMasker_Tools/RepeatMasker/RepeatMasker
LTR=/home/why/Tools/RepeatMasker_Tools/LTR_Finder/source/ltr_finder
LTR_P=/home/why/Tools/RepeatMasker_Tools/LTR_FINDER_parallel/LTR_FINDER_parallel

genome=/home/why/maker_annotation/RA_repeatmasking/New_polish/RA.final_polished.fasta
species=RA

$LTR_P -seq $genome -threads 20 -harvest_out #> $specoes.LTRFINDER.scn

LTR_RETRIEVER=/home/why/Tools/RepeatMasker_Tools/LTR_retriever/LTR_retriever
repeatmasker="/home/why/Tools/RepeatMasker_Tools/RepeatMasker/"
cdhit="/home/why/Tools/RepeatMasker_Tools/cdhit"

genome=/home/why/maker_annotation/RA_repeatmasking/New_polish/RA.final_polished.fasta
scn=RA.final_polished.fasta.finder.combine.scn


$LTR_RETRIEVER -genome $genome -inharvest $scn -threads 20 -repeatmasker $repeatmasker -cdhit_path $cdhit

#Final (Remove_oldmasked_then_remasked)
cat ../De_Novo/RM_587686.ThuAug111114132022/consensi.fa.classified ../RepBase_Dfam/Percomorphaceae_TE_lib.fa ../PM_TE/Periophthalmus_magnuspinnatus_transposable_elements ../LTR/RA.final_polished.fasta.mod.LTRlib.fa > RA_all_four_TE_lib.fa
/home/why/Tools/RepeatMasker_Tools/RepeatMasker/RepeatMasker -s -pa 10 -gff -engine rmblast -lib RA_all_four_TE_lib.fa -xsmall RA.final_polished_unmasked.sorted.fa

#Final2 (Remove_oldmasked_then_remasked)
cat ../De_Novo/RM_587686.ThuAug111114132022/consensi.fa.classified ../RepBase_Dfam/Percomorphaceae_TE_lib.fa ../PM_TE/Periophthalmus_magnuspinnatus_transposable_elements > RA_three_TE_lib.fa
/home/why/Tools/RepeatMasker_Tools/RepeatMasker/RepeatMasker -s -pa 10 -gff -engine rmblast -lib RA_three_TE_lib.fa -xsmall RA.final_polished_unmasked.sorted.fa