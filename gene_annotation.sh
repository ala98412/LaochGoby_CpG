braker.pl --useexisting --species=Rhyacichthys_aspro_RNA --genome=RA.genome.fa\
 --bam=gill.sorted.bam,brain.sorted.bam,eyeB.sorted.bam,liver1.sorted.bam,nose.sorted.bam,muscle.sorted.bam --softmasking --workingdir=RA_RNA_out\
 --cores 10 -gff3 --UTR=on --nocleanup
braker.pl --species=Rhyacichthys_aspro_RNA_2 --genome=RA.genome.fa --bam=gill.sorted.bam,brain.sorted.bam,eyeB.sorted.bam,liver1.sorted.bam,nose.sorted.bam,muscle.sorted.bam --softmasking --workingdir=RA_RNA_out --cores 15 -gff3 --UTR=on --nocleanup
/home/why/BREAKER_annotation/TSEBRA/bin/fix_gtf_ids.py --gtf RA_RNA_out/braker.gtf --out RA_RNA_out_fixed.gtf

cat BP_PM_SO.pep.fa Betta_splendens.pep.fa Scophthalmus_maximus.pep.fa Oryzias_latipes.pep.fa Nothobranchius_furzeri.pep.fa Gasterosteus_aculeatus.pep.fa Larimichthys_crocea.pep.fa Takifugu_rubripes.pep.fa Syngnathus_acus.pep.fa Thalassophryne_amazonica.pep.fa > ../Percomorphaceae.new.fa
/home/why/BREAKER_annotation/ProtHint/bin/prothint.py --threads 10 RA.genome.fasta Percomorphaceae.new.fa
braker.pl --useexisting --species=Rhyacichthys_aspro_pro --genome=RA.genome.fa\
 --hints=RA_ProHint/prothint_augustus.gff --softmasking --workingdir=RA_pro_out\
 --cores 10 -gff3 --nocleanup
/home/why/BREAKER_annotation/TSEBRA/bin/fix_gtf_ids.py --gtf RA_pro_out/braker.gtf --out RA_pro_out_fixed.gtf

/home/why/BREAKER_annotation/TSEBRA/bin/tsebra.py -g RA_RNA_out_fixed.gtf,RA_pro_out_fixed.gtf -c ../TSEBRA/config/pref_braker1.cfg -e RA_RNA_out/hintsfile.gff,RA_pro_out/hintsfile.gff -o RA_combined.gtf
/home/why/BREAKER_annotation/TSEBRA/bin/rename_gtf.py --gtf RA_combined.gtf --out RA_combined_new.gtf #important
/home/why/BREAKER_annotation/Augustus/scripts/gtf2gff.pl < RA_combined_new.gtf --out=RA_combined.gff3 --gff3
/home/why/BREAKER_annotation/Augustus/scripts/getAnnoFastaFromJoingenes.py -g RA.genome.fa -f RA_combined_new.gtf -o "TSEBRA_PrefBraker1_RA"