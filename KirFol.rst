# reStructuredText

## SNP filtering; kirfol/vcfs/kf/indv_vcfs
### gatk.ScoreCNN.qsub.script  # scorerealign.vcf
	-grep -w "alignment_artifact" > realign.vcf.aln.mask
	-bcftools view -f "." > realign.flt1.vcf
	-CNNScoreVariants < realign.flt1.vcf > realign.flt1.ScoreCNN.vcf
	-FilterVariantTranche < realign.flt1.ScoreCNN.vcf > realign.flt1.ScoreCNN.Tranche.vcf
	-sort -m -k2,2n <(gunzip -c realign.flt1.ScoreCNN.Tranche.vcf.gz) <(gunzip -c realign.vcf.aln.mask.gz) | gzip -c > realign.flt1.ScoreCNN.Tranche.aln.vcf.gz
		>rm *realign.vcf.aln.mask

### population VCFs; kirfol/vcfs/kf
	-KirFol.AfunF3.Contig0.3R.pop.vcf
	-KirFol.AfunF3.Contig0.3L.pop.vcf
	-KirFol.AfunF3.Contig1.2R.pop.vcf
	-KirFol.AfunF3.Contig1.2L.pop.vcf
	-KirFol.AfunF3.Contig2.X.pop.vcf

### Create masking file
	-vcf2mask.py -f *realign.flt1.ScoreCNN.Tranche.aln.vcf.gz -o KirFol.VCFmask.$CHR.txt
	-applymask2vcf.py -f KirFol.AfunF3.$CONTIG.$CHR.pop.vcf.gz -m KirFol.mask.txt.gz -o KirFol.AfunF3.Contig{0,1,2}.${CHR}.pop.flt.vcf

### post-processing; remove indels, gtDP gtGQ
	-vcftools --gzvcf KirFol.AfunF3.Contig{0,1,2}.$CHR.pop.flt.vcf.gz --remove-indels --minDP 10 --minGQ 30 --recode --stdout | gzip -c > KirFol.AfunF3.Contig{0,1,2}.${CHR}.pop.flt.HQ.snps.vcf.gz

	X: kept 14459274 out of a possible 16211358 Sites
	2R:
	2L:
	3R:
	3L:

### missing mask for fastas
	-missmask.py -f KirFol.AfunF3.Contig{0,1,2}.${CHR}.pop.flt.HQ.snps.vcf.gz
		>$vcfFile.genome.mask.bed   # use for masking base-fastafile, these are missing from all AFTER individual filtering
		>Individual.mask.$CHR.txt.gz  # use for individual masks, sites where not all missing but individuals have missing gt
		!assume no more filtering on individuals!
	
	X: 
	2R:
	2L:
	3R:
	3L:

### nonaccessible regions + total callable sites
	-vcftools --gzvcf KirFol.AfunF3.Contig{0,1,2}.$CHR.pop.flt.HQ.snps.vcf.gz --exclude-bed nonaccessible.merge.bed --recode --stout | gzip -c > KirFol.AfunF3.Contig{0,1,2}.$CHR.pop.flt.HQ.snps.mask.vcf.gz
		>total sites with at least 1 called individual (Het + HomR + HomA)

### remove samples with > X% missing data
	-vcf_to_zarr : scikit-allel
		?individual missingness?
		?how many sites w/ X% missing?
		?how many individuals to remove to maximize site coverage?
		!total sites should not have changed!
		>pop.flt.HQ.snps.mask.143.vcf.gz  # add number individuals
	-vcftools --gzvcf --site-missingness  # just some stats here, prob same from sk-allel

### total callable sites with no missing	
	vcftools --gzvcf KirFol.AfunF3.Contig{0,1,2}.$CHR.pop.flt.HQ.snps.mask.143.vcf.gz --max-missing 1 --out totalCallable.Nomissing
		!this is total number of callable sites with no missing data!

### phasing
	-shapeit4 w/recomb.map from LDJump, PS option
	-addphase2vcf.py, file will have phased and non-phased
	?what percent missing data is allowed?
	?tri-nucleotide sites?  # code individuals as missing?, or depends on percentage?
	?which sites are removed before phasing?  # need to keep these to mask
	!Sites that cant be phased need to be masked in phased fasta!

### phased fasta
	-vcf2fasta.py 
		!reference fasta masked with nonaccessible + Genome.mask.bed
		!individual fasta masked with Individual.mask.$CHR.txt.gz
		!add phased mask, sites that are called but not phased!

KEEP
*pop.flt.vcf.gz  # has indels, no filtering
>*pop.flt.HQ.snps.vcf.gz  # has snps, GQ and DP filtering
>*pop.flt.HQ.snps.mask.vcf.gz  # nonaccessible mask
>*pop.flt.HQ.snps.mask.143.vcf.gz  # remove individuals with too much missing data
*pop.flt.HQ.snps.mask.143.var.vcf.gz  # variable sites only, no HomR, non-ref-ac 1, all are "/"
*pop.flt.HQ.snps.mask.143.var.mix.vcf.gz  # variable sites only, breaks in synteny but "|" added, so mix of "|" and "/" if 0/0 or 1/1 can be phased for that individual
*pop.flt.HQ.snps.mask.143.var.phased.vcf.gz  # phased sites only


## ancestral state estimations
vcfmerge subgroup & kirfol
ancestral state using keightly
	>FOO.miss.anc.vcf

## analysis
allel - LD thinning


export GATK_LOCAL_JAR=/afs/crc.nd.edu/user/s/ssmall2/programs_that_work/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar
source gatk
export PATH=$HOME/anaconda2/bin:$PATH
export GATK_SPARK_JAR=/afs/crc.nd.edu/user/s/ssmall2/programs_that_work/gatk-4.1.0.0/gatk-package-4.1.0.0-spark.jar