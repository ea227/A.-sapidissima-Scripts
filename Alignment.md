Quality controlled fastq files were  aligned to the reference using bwa mem:

```
bwa mem \ #Aligning libraries to reference genome 
   -M \ 
   -R "@RG\tID:${out}\tSM:${out}\tPL:illumina\tLB:run1" \ 
   -t 16 \ 
   $ref \ 
   $r1 \ 
   $r2 | \ 
```
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
