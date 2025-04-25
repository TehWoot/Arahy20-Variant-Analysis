# Introduction

A variant analysis was conducted on QTL *qPSIIB10*. The *qPSIIB10* region has been shown as a locus possibly conferring aflatoxin resistance in peanut and is present across multiple environments. In the present study, susceptible and resistant varieties, Florida-07 and ICG1471, were mapped to cultivated peanut, Tifrunner, to identify variants and candidate genes possibly contributing to aflatoxin resistance.

# Available Data

The data available in this repository is listed as follows:

"Flor07_chr20_snp_table.txt" & "ICG1471_chr20_snp_table.txt"
> Variants called using the SnpEff annotation tool on chromosome 20.

"Flor07_SnpEff_Chr20FilteredQual30.html" & "ICG1471_SnpEff_Chr20FilteredQual30.html"
> SnpEff weblink data of variants present on chromosome 20.

"Flor07_snpEff_summary.html" & "ICG1471_snpEff_summary.html"
> SnpEff weblink data of variants present in the Tifrunner genome.

"updatedduip.fa.txt2.fa_mq10.sam.txt3.bed"
> Differentially expressed genes in response to *Aspergillus flavus* infection mapped by Walid Korani.

The data available in [this Google Drive](https://drive.google.com/drive/folders/1KtHvaQtU1dmY667oeUHx70eUEflsZOD3?usp=sharing) is listed as follows.

"Flor07align_srt.bam" & "ICG1471align_srt.bam"
> Aligned sequences of Florida-07 and ICG1471.

"Flor07calls.vcf" & "ICG1471calls.vcf"
> Variant calls of Florida-07 and ICG1471.

"Flor07calls.bcf" & ICG1471calls.bcf"
> Variant calls of Florida-07 and ICG1471.





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
Variants were viewed using IGV 2.18.4.

# SnpEff Analysis

A SnpEff analysis was conducted on the mapped variants from Florida-07 and ICG1471 on chromosome 20. SnpEff is a genetic variant annotation and functional effect prediction toolbox (Cingolani et al., 2012). Both the Florida-07 and ICG1471 cultivars were mapped and had variants called to allow for the detection of SNPs in this analysis. Building a SnpEff database and filtering was done according to [this tutorial](https://www.youtube.com/watch?v=-rmreyRAbkE&ab_channel=RodrigoBaptista). Files containing variants ranked as high, moderate, low, and modifier impact, as well as detailed SnpEff web links,  are listed in the repository.

```
# SnpEff Annotation with Filtered .vcfs

java -jar /snpeff_2/snpEff.jar ara_hypo vcfs/ICG1471calls_chr20_filtered.vcf > ICG1471calls_chr20_filtered.ann.vcf
```

```
# Filtering Annotated SnpEff.vcf

more ICG1471calls_chr20_filtered.ann.vcf | grep -v "#" | grep "PASS" | cut -f 1,2,8 | sed 's/|/\t/g' | cut -f 1,2,4,5,7 | grep -v "intergenic\|stream\|UTR" > ICG1471_snp_table
```
The SnpEff data collected using variants called in the peanut genome supports previous work done by XYZXYZXYZ . Further showing that the susceptible variety peanut, Florida-07, is closely related to 


<img width="444" alt="ICG1471SnpEffChart" src="https://github.com/user-attachments/assets/be294ccc-0446-4af8-8d87-ed014d9abc2c" />


<img width="371" alt="Flor07SnpEffChart" src="https://github.com/user-attachments/assets/c8063b58-324e-499f-bdf1-1db9162eca31" />


# Differentially Expressed Genes Mapped to Tifrunner

Korani et al. (2018) identified differentially expressed genes from the susceptible and resistant peanut cultivars Florida-07 and ICG1471. These DE genes were mapped to TifrunnerV2. The mapped DE genes were provided by Walid Korani. The differentially expressed genes were viewed in IGV 2.18.4.

<img width="1920" alt="DEGenesOutsideqPSBIIB10" src="https://github.com/user-attachments/assets/e487eb77-310b-4878-8db0-4fc53c329208" />

# Tools Used









