# Introduction

A variant analysis was conducted on QTL *qPSIIB10*. The *qPSIIB10* region has been shown as a locus possibly conferring aflatoxin resistance in peanut and is present across multiple environments. In the present study, resistant and susceptible varieties, ICG1471 and Florida-07, were mapped to cultivated peanut, Tifrunner, to identify variants and candidate genes possibly contributing to aflatoxin resistance.

# Mapping and Variant Calling

The HiFi sequences of Florida-07 and ICG1471 were mapped to TifrunnerV2 using minimap2 2.26. Variants were called using samtools 1.16.1 and bcftools 1.18.

```
minimap2 -a arahy.Tifrunner.gnm2.J5K5.genome_main.fna ICG1471_HIFI.fastq > ICG1471align.sam

samtools view -bS -@ 8 ICG1471align.sam > ICG1471align.bam
samtools sort -@ 8 ICG1471align.bam -o ICG1471align_srt.bam
samtools index -c -@ 8 ICG1471align_srt.bam
```

```
bcftools mpileup -Ou -f arahy.Tifrunner.gnm2.J5K5.genome_main.fna ICG1471align_srt.bam | bcftools call -mv -Ob -o ICG1471calls.bcf

bcftools convert -Oz -o ICG1471calls.vcf ICG1471calls.bcf
```

```
# Filtering Variants
bcftools filter -i 'QUAL>=30' ICG1471calls.vcf | grep -v -c '^#'
```
Variants were viewed using IGV 2.18.4

# SnpEff Analysis

A SnpEff analysis was conducted on the mapped variants from Florida-07 and ICG1471 on chromosome 20. SnpEff is a genetic variant annotation and functional effect prediction toolbox (Cingolani et al., 2012). Both the Florida-07 and ICG1471 cultivars were mapped and had variants called to allow for the detection of SNPs in this analysis. Files containing variants ranked as high, moderate, low, and modifier impact are listed in the repository.

```
# SnpEff Annotation with Filtered .vcfs

java -jar /snpeff_2/snpEff.jar ara_hypo vcfs/ICG1471calls_chr20_filtered.vcf > ICG1471calls_chr20_filtered.ann.vcf
```

```
# Filtering Annotated SnpEff.vcf

more ICG1471calls_chr20_filtered.ann.vcf | grep -v "#" | grep "PASS" | cut -f 1,2,8 | sed 's/|/\t/g' | cut -f 1,2,4,5,7 | grep -v "intergenic\|stream\|UTR" > ICG1471_snp_table
```








