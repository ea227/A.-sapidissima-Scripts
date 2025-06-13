# **Alignment and Variant Calling**
### Step 1: Map reads
All cleaned, paired-end reads were aligned to the reference genome using bwa mem. This was repeated separately for all three libraries as one large pipeline.

_mapping with bwa mem_
```bash
# Setup variables; repeated for each SRR library
SRR=SRR7973879
threads=48 # CPU threads to use
rg="@RG\tID:${SRR}\tSM:${SRR}\tPL:illumina\tLB:run1"

# Run bwa mem
bwa mem \
      -M \
      -p \
      -R "$rg" \
      -t ${threads} \
      Asap.fa \
      ${SRR}.finalclean_1.fq.gz \
      ${SRR}.finalclean_2.fq.gz | the output was piped into samtools
```

- _Parameters Explained:_
  - ***-M*** :: mark shorter split hits as secondary
  - ***-p*** :: smart pairing (ignoring in2.fq), in other words, data input are interleaved fastq
  - ***-R*** :: Read group ID line (i.e. sample name)
  - ***-t*** :: number of cpus to use

_initial processing with samtools_
```bash
# Process with samtools
samtools sort \
   -T ${SRR}.tmp \
   -n \
   -@ ${threads} \
   - | \
   samtools fixmate \
      -@ ${threads} \
      -m \
      - - | \
      samtools sort \
         -T ${SRR}.tmp2 \
         -O bam \
         -@ ${threads} \
         - | \
         samtools markdup \
             -T ${SRR}.tmp3 \
             -O bam \
             -@ ${threads} \
             - ${SRR}.sorted.bam
```

- _Parameters Explained:_
  - ***sort -n*** :: sort BAM file numerically
  - ***fixmate -m*** :: fixmates and add mate score tag
  - ***markdup*** :: mark PCR/optical duplicates for later removal.

_post-process BAM file_
```bash
# Process BAM and get stats
samtools index -@ ${threads} ${SRR}.sorted.bam
samtools stats -@ ${threads} ${SRR}.sorted.bam > ${SRR}.bamstats
```

### Step: Merge BAM files

The resulting three BAM files were merged with `samtools merge` into a single BAM file:

```bash
threads=48 # CPU threads to use

# Merge bamfiles
samtools merge \
      --threads ${threads} \
      -o shad.merged.bam \
      SRR7973879.sorted.bam \
      SRR7973880.sorted.bam \
      SRR7973881.sorted.bam

# Index
samtools index shad.merged.bam
samtools stats -@ ${threads} shad.merged.bam > shad.merged.bamstats
```


# XXXXXXXXXXXXXXXXXXXXX

```
Variants were called with bcftools v1.10.2 mpileup, filtered, and output as a gzipped VCF:
```
bcftools \  
        mpileup \ 
        -C 50 \ 
        -q 20 \ 
        -Q 20 \ 
        --threads 24 \ 
        -Ou \ 
        --ignore-RG \ 
        --per-sample-mF \ 
  --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \ 
        -f /${path} \ 
        shad.merged.bam | \ 
        bcftools \ 
           call \ 
           --threads 24 \ 
           --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \ 
           -mv -Ou | \ 
        bcftools \ 
           filter \ 
           --threads 24 \ 
           -sFAIL -e'QUAL < 30 || INFO/MQ <= 30 || INFO/MQBZ < -9 || INFO/RPBZ < -5 || INFO/RPBZ > 5 || INFO/BQBZ < -5 || INFO/DP4[3]+INFO/DP4[4] <= 2' \ 
           -g10 \ 
           -G10 \ 
           -Ov | \ 
        bgzip -c > shad.vcf.gz 
```
bcftools was then used to index the VCF, generate stats, and filter the VCF further:
```
bcftools index --threads 24 shad.vcf.gz 
bcftools stats shad.vcf.gz > shad.SNPs.stats 
bcftools view -f 'PASS' -O v shad.vcf.gz | bgzip -c > shad.filtered.vcf.gz 
bcftools stats -f "PASS" shad.filtered.vcf.gz > shad.filtered.SNPs.stats
```
Variant information was done with SnpEff:
```
java -jar /home/rfitak/PROGRAMS/snpEff/snpEff.jar ann \
        -c /home/rfitak/PROGRAMS/snpEff/snpEff.config \
        -v aloSap \
        shad.filtered.vcf.gz | \
        bgzip -c > shad.filtered.vann.vcf.gzip
```

Variant information was accessed for 20000 bp windows using vcftools v0.1.17, and output as a .tsv file: 
```
vcftools --gzvcf shad.vcf.gz --SNPDensity 20000 --out filename
```

