# Introduction

A variant analysis was conducted on QTL *qPSIIB10*. The *qPSIIB10* region has been shown as a locus possibly conferring aflatoxin resistance in peanut and is present across multiple environments. In the present study, resistant and susceptible varieties, ICG1471 and Florida-07, were mapped to cultivated peanut, Tifrunner, to identify variants and candidate genes possibly contributing to aflatoxin resistance.

# Mapping and Variant Calling

The HiFi sequences of Florida-07 and ICG1471 were mapped to TifrunnerV2 using minimap2 2.26. Variants were called using samtools 1.16.1 and bcftools 1.18.
Variants were viewed using IGV 2.18.4

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



