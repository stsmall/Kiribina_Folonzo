# KirFol

* Author: S.T. Small

Overview of SNP calling and filtering pipeline for An. funestus

## Inputs

* Fasta Ref
* PE illumina reads
* mapping for readgroups

## Outputs

* Mapped BAM File w/ readgroup
* basic statistics from QualiMap
* gVCF

## Software
* bwa v0.7.17-r1188
* samtools v1.3.1
* picard v2.18.16
* qualimap v2
* GATK v4.1.2
* bedtools v2.28.0

## Pipeline

### Align bwa mem

```
bwa mem -q -M -R "readgroup_string" ref read1 read2 | samtools view -bS - | samtools sort -@ 5 - > output.bam
```

Notes: 

* -q "don't modify mapQ of supplementary alignments"
* -M mark shorter split hits as secondary

### MarkDuplicates

```
picard.jar MarkDuplicates
```

### Validate

```
picard.jar ValidateSamFile
```

### Index

```
samtools index file.bam
```

### Statistics

```
samtools flagstat file.bam
```

```
genomeCoverageBed -d -ibam file.bam
```

```
qualimap bamqc -bam file.bam -outdir
qualimap_results -outformat pdf
```



### SNP calling & Filtering

#### create truth set by hard filtering

> Keep called by both strelka and HaplotypeCaller

1. ```strelka2 [option]```  
2. ```gatk HaplotypeCaller -R ref -I Sample.bam --emit-ref-confidence GVCF --heterozygosity 0.01 --indel-heterozygosity 0.001 --min-base-quality-score 17 --bam-output Sample.HC.bam -L contig -O Sample.g.vcf```  
3. ```gatk GenotypeGVCFs -R AfunF3.primary.fasta --heterozygosity 0.01 --indel-heterozygosity 0.001 -new-qual -V g.vcf -O output.vcf```
4. ```gatk FilterAlignmentArtifacts -R ref -V vcf -I HC.bam --bwa-mem-index-image ref.fasta.img -O realign.vcf```
5. ```gatk SplitVcfs -I realign.vcf --SNP_OUTPUT SNP.vcf --INDEL_OUTPUT INDEL.vcf --STRICT false```  
6. ```gatk VariantFiltration -R ref -V SNP.vcf -O SNP.flt.vcf --set-filtered-genotype-to-no-call true --filter-name "LowQual" --filter-expression "QD < 2.0 || MQ0 < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || FS > 60.0 || SOR > 3.0"```  
7. ```gatk VariantFiltration -R ref -V INDEL.vcf -O INDEL.flt.vcf --set-filtered-genotype-to-no-call true --filter-name "LowQual" --filter-expression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0 || SOR > 10.0"```
8. ```gatk MergeVcfs -I SNP.flt.vcf -I INDEL.flt.vcf -O realign.flt.vcf```


#### Train CNN for SNP filtering



1. ```gatk CNNVariantWriteTensors -R ref -V realign.flt.vcf -truth-vcf truth.TrainCNN.vcf -truth-bed accessible.bed -tensor-type read_tensor -bam-file HC.bam -output-tensor-dir ${contig}-tensor```
2. ```gatk CNNVariantTrain -input-tensor-dir ${contig}-tensor/ -tensor-type read_tensor -model-name ${contig}_2d_model```

**Note**

- accessible.bed is callable sites from SNPable and using GFF to mask remove simple and tandem repeats  
- GFF is available on VectorBase for AfunF3 assembly


#### Population Genotyping


> The problem is always back-filling from IND to POP. FilterAlignArt runs much faster on just SNPs and the CNN is trained only on SNPs.
> I am not sure what would happen if one used all sites, other than it would take a long time. A possibility might be to split the vcf
> containing all the sites to SNPs and HomR and then do the scoring. I trained and filtered on a SNP only vcf. I then regenotyped all as
> a population using filtered gvcfs. This was the long way around, but I was testing the CNN and AlignArt at that time.
> Regardless, I think that the below option is more flexible since you keep the VCF in the case of filtering changes.


1. ```gatk HaplotypeCaller -R ref -I Sample.bam --emit-ref-confidence GVCF --heterozygosity 0.01 --indel-heterozygosity 0.001 --min-base-quality-score 17 --bam-output Sample.HC.bam -L contig -O Sample.g.vcf```  
2. ```gatk GenotypeGVCFs -R ref --heterozygosity 0.01 --indel-heterozygosity 0.001 -new-qual -all-sites -V Sample.g.vcf -O Sample.all.vcf```  
3. ```SPLIT Variants > [Sample.SNP.vcf, Sample.HomR.vcf]```  
4. ```gatk FilterAlignmentArtifacts -R ref -V Sample.SNP.vcf -I HC.bam --bwa-mem-index-image ref.fasta.img -O realign.vcf```  
5. ```gatk CNNScoreVariants -R ref -I HC.bam -V realign.vcf -O ScoreCNN.vcf -tensor-type read_tensor -architecture 2d_model.json -weights 2d_model.hd5 --inter-op-threads 5 --intra-op-threads 5```  
6. ```gatk FilterVariantTranches -V ScoreCNN.vcf --resource truth.ScoreCNN.vcf --info-key CNN_2D -snp-tranche 90.0 -snp-tranche 95.0 -snp-tranche 97.0 -snp-tranche 99.0 -snp-tranche 99.90 -indel-tranche 90.0 -indel-tranche 95.0 -indel-tranche 99.90 --invalidate-previous-filters true -O tranche.vcf```  
7. ```REMERGE tranche.vcf Sample.HomR.vcf > Sample.all.flt.vcf```
8. ```gatk VariantFiltration [options]```  
**To get PopSet** 
9. ```gatk CombineVariants OR picard.jar MergeVcfs```  
10. ```python getVCFStats.py```  




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

### post-processing; remove indels, gtDP gtGQ ? module load gcc/6.2.0
	-vcftools --gzvcf KirFol.AfunF3.Contig{0,1,2}.$CHR.pop.flt.vcf.gz --remove-indels --minDP 10 --minGQ 30 --recode --stdout | gzip -c > KirFol.AfunF3.Contig{0,1,2}.${CHR}.pop.flt.HQ.snps.vcf.gz

	X : kept 14459274 out of a possible 16211358 Sites
	2R: kept 46115694 out of a possible 50348497 Sites
	2L: kept 36683357 out of a possible 40708819 Sites
	3R: kept 42394418 out of a possible 45933782 Sites
	3L: kept 36524906 out of a possible 39979476 Sites
			176177649				   193181932

### missing mask for fastas
	-missmask.py -f KirFol.AfunF3.Contig{0,1,2}.${CHR}.pop.flt.HQ.snps.vcf.gz
		>$vcfFile.genome.mask.bed   # use for masking base-fastafile, these are missing from all AFTER individual filtering
		>Individual.mask.$CHR.txt.gz  # use for individual masks, sites where not all missing but individuals have missing gt
		!assume no more filtering on individuals!
	
### nonaccessible regions + total callable sites
	-vcftools --gzvcf KirFol.AfunF3.Contig{0,1,2}.$CHR.pop.flt.HQ.snps.vcf.gz --remove-indv La_02-07553_Folonzo --max-missing-count 202 --exclude-bed /scratch365/ssmall2/reference_fasta_index/AnfunSG_refs/AfunF3/nonaccessible.merge.bed --recode --stdout | gzip -c > KirFol.AfunF3.Contig{0,1,2}.$CHR.pop.flt.HQ.snps.mask.vcf.gz
		#remove sites with all missing
		>total sites with at least 1 called individual (Het + HomR + HomA)

	X : kept 11202630 out of a possible 14459274 Sites
	2R: kept 35514376 out of a possible 46115694 Sites
	2L: kept 28478106 out of a possible 36683357 Sites
	3R: kept 28946333 out of a possible 42394418 Sites
	3L: kept 20205402 out of a possible 36524906 Sites
            124346847				   176177649

### remove samples with > X% missing data  *** TODO
	vcftools --gzvcf --missing-indv
	vcftools --gzvcf --missing-site
	-vcf_to_zarr : scikit-allel; add plots
		?how many sites w/ X% missing?
		?how many individuals to remove to maximize site coverage?
		!total sites should not have changed!
		>pop.flt.HQ.snps.mask.143.vcf.gz  # add number individuals

### total callable sites with no missing	
	vcftools --gzvcf KirFol.AfunF3.Contig{0,1,2}.$CHR.pop.flt.HQ.snps.mask.143.vcf.gz --max-missing 1 --out totalCallable.Nomissing
		!this is total number of callable sites with no missing data!
	X : ?
	2R: ?
	2L: ?
	3R: ?
	3L: ?

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

NOTES
*pop.flt.vcf.gz  # has indels, no filtering
>*pop.flt.HQ.snps.vcf.gz  # has only snps, GQ and DP filtering
>*pop.flt.HQ.snps.mask.vcf.gz  # nonaccessible mask
>*pop.flt.HQ.snps.mask.143.vcf.gz  # remove individuals with too much missing data
*pop.flt.HQ.snps.mask.143.var.vcf.gz  # variable sites only, no HomR, non-ref-ac 1, all are "/"
*pop.flt.HQ.snps.mask.143.var.mix.vcf.gz  # variable sites only, breaks in synteny but "|" added, so mix of "|" and "/" if 0/0 or 1/1 can be phased for that individual
*pop.flt.HQ.snps.mask.143.var.phased.vcf.gz  # phased sites only


### ancestral state estimations
vcfmerge subgroup & kirfol
ancestral state using keightly
	>FOO.miss.anc.vcf add as AA

### analysis

## population structure
allel - LD thinning

### jupyter notebook
~/anaconda3/bin/jupyter-notebook --no-browser --port=8777
ssh -f ssmall2@rosalind.crc.nd.edu -L 8777:localhost:8777 -N
paste link in browser

KEEP
*realign.flt1.ScoreCNN.Tranche.aln.vcf.gz  # move to backup drive
*pop.flt.vcf.gz  # move to backup

*pop.flt.HQ.snps.mask.143.var.vcf.gz
*pop.flt.HQ.snps.mask.143.var.mix.vcf.gz
*pop.flt.HQ.snps.mask.143.var.phased.vcf.gz

DELETE
rm -f *realign.vcf.gz
rm -f *realign.flt1.ScoreCNN.Tranche.vcf.gz
rm -f *realign.vcf.aln.mask.gz
rm -rf indiv_gvcfs
rm -f *pop.flt.HQ.snps.vcf.gz
rm -f *pop.flt.HQ.snps.mask.vcf.gz
rm -f *pop.flt.HQ.snps.mask.143.vcf.gz
