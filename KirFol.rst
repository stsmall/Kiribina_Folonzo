# reStructuredText

## kirfol/vcfs/kf/indv_vcfs
gatk.ScoreCNN.qsub.script
realign.vcf
	#  grep -w "alignment_artifact" > realign.vcf.aln.mask
	#  bcftools view -f "." > realign.flt1.vcf
CNNScoreVariants < realign.flt1.vcf > realign.flt1.ScoreCNN.vcf
FilterVariantTranche < realign.flt1.ScoreCNN.vcf > realign.flt1.ScoreCNN.Tranche.vcf
realign.flt1.ScoreCNN.Tranche.alnart.vcf < realign.vcf.aln.mask < realign.flt1.ScoreCNN.Tranche.vcf


## kirfol/vcfs/kf
KirFol.AfunF3.Contig0.3R.pop.vcf
KirFol.AfunF3.Contig0.3L.pop.vcf
KirFol.AfunF3.Contig1.2R.pop.vcf
KirFol.AfunF3.Contig1.2L.pop.vcf
KirFol.AfunF3.Contig2.X.pop.vcf

## Masking File
vcf2mask.py -f *realign.flt1.ScoreCNN.Tranche.alnart.vcf -o KirFol.VCFmask.txt
applymask2vcf.py -f *.pop.vcf -m KirFol.mask.txt -o KirFol.AfunF3.Contig{0,1,2}.${CHR}.pop.flt.vcf
missmask.py -f KirFol.AfunF3.Contig{0,1,2}.${CHR}.pop.flt.vcf
	>vcfFile.genome.mask.bed
	>Individual.mask.$CHR.txt.gz

## post-processing
vcftools --gzvcf KirFol.AfunF3.Contig{0,1,2}.$CHR.pop.flt.vcf -exlclude nonaccessible.bed --non-ref-ac 1 --recode -o FOO
python remove_missinds.py
ancestral state using keightly

## phasing
shapeit2 w/ recomb.map from LDJump fun and pir gatk.bam
python addphase2vcf.py

## analysis - formats
Zarr to scikit-allel
plink
vcf


