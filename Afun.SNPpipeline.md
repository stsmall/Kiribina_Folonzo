# An. funestus SNP Calling

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