## **Alignment and Variant Calling**

All cleaned, paired-end reads were aligned to the reference genome using bwa mem.

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
      ${SRR}_1.... \
      ${SRR}_2.... | the output was piped into samtools
```

- _Parameters Explained:_
  - ***-M*** :: mark shorter split hits as secondary
  - ***-p*** :: smart pairing (ignoring in2.fq), in other words, data input are interleaved fastq
  - ***-t*** :: number of cpus to use

_initial processing with samtools_
```bash
# Process with samtools
   # ${threads} = # CPU threads to use
   # ${out} = output file basename
samtools sort \
   -T ${out}.tmp \
   -n \
   -@ ${threads} \
   - | \
   samtools fixmate \
      -@ ${threads} \
      -m \
      - - | \
      samtools sort \
         -T ${out}.tmp2 \
         -O bam \
         -@ ${threads} \
         - | \
         samtools markdup \
             -T ${out}.tmp3 \
             -O bam \
             -@ ${threads} \
             - ${out}.sorted.bam
```

- _Parameters Explained:_
  - ***sort -n*** :: sort BAM file numerically
  - ***fixmate -m*** :: fixmates and add mate score tag
  - ***markdup*** :: mark PCR/optical duplicates for later removal.

_cleaning the BAM file_
```bash
# Clean up the bam file
samtools view \
      -b \
      -h \
      -q 20 \
      -f 0x2 \
      -F 0x4 \
      -F 0x8 \
      -F 0x400 \
      -@ ${threads} \
      -o ${out}.clean.sorted.bam \
      ${out}.sorted.bam
```

- _Parameters Explained:_
  - ***-b*** :: output BAM format
  - ***-h*** :: Include header in SAM output
  - ***-q*** :: remove reads with mapping quality < 20
  - ***-f 0x2*** :: keep reads mapped in proper pair
  - ***-F 0x4*** :: remove unmapped reads
  - ***-F 0x8*** :: remove read when its mate is unmapped
  - ***-F 0x400*** :: remove read when mate is mapped to the reverse strand, or not the primary alignment.

_post-process BAM file_
```bash
# Process BAM and get stats
   # ${threads} = # CPU threads to use
   # ${out} = output file basename
samtools index -@ ${threads} ${out}.clean.sorted.bam
samtools stats -@ ${threads} ${out}.clean.sorted.bam > ${out}.bamstats
```





Quality controlled fastq files were  aligned to the reference using the Burrows-Wheeler Aligner BWA MEM:

```
bwa mem \ #Aligning libraries to reference genome 
   -M \ 
   -R "@RG\tID:${out}\tSM:${out}\tPL:illumina\tLB:run1" \ 
   -t 16 \ 
   $ref \ 
   $r1 \ 
   $r2 | \ 
```
BWA MEM parameters:

-M: Mark shorter split hits as secondary

-R: Read group ID line (i.e. sample name)

-t: number of threads

samtools v 1.10.2 was used to sort BAM files numerically, fixmates to add mate score tag and markdup to identify  duplicates:
```
samtools sort \ 
   -T ${out}.tmp \ 
   -n \ 
   -@ 16 \ 
   - | \
   samtools fixmate \ 
      -@ 16 \ 
      -m \ 
      - - | \ 
      samtools sort \ 
         -T ${out}.tmp2 \ 
         -O bam \ 
         -@ 16 \ 
         - | \
         samtools markdup \ 
             -T ${out}.tmp3 \ 
             -O bam \ 
             -@ 16 \ 
             - ${out}.sorted.bam
```
The resulting BAM files were merged with samtools merge:
```
samtools merge \ 
        --threads 24 \ 
        -o shad.merged.bam \ 
        Path$(SRR)  \ 
        Path$(SRR)   \     
        Path${SRR}
```
## **Processing**

The merged alignment was then processed:

```
samtools index shad.merged.bam #Indexing merged bam file 
samtools stats -@ 24 shad.merged.bam > shad.merged.bamstats  
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

