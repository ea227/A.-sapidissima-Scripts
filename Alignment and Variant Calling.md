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
      -R "$rg" \
      -t ${threads} \
      Asap.fa \
      ${SRR}.finalclean_1.fq.gz \
      ${SRR}.finalclean_2.fq.gz | the output was piped into samtools
```

_Parameters Explained:_
  - ***-M*** :: mark shorter split hits as secondary
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

_Parameters Explained:_
  - ***sort -n*** :: sort BAM file numerically
  - ***fixmate -m*** :: fixmates and add mate score tag
  - ***markdup*** :: mark PCR/optical duplicates for later removal.

_post-process BAM file_
```bash
# Process BAM and get stats
samtools index -@ ${threads} ${SRR}.sorted.bam
samtools stats -@ ${threads} ${SRR}.sorted.bam > ${SRR}.bamstats
```

The various summary statistics of the mapping for each set of reads can be found here:
- [SRR7973879](./data/SRR7973879.bamstats)
- [SRR7973880](./data/SRR7973880.bamstats)
- [SRR7973881](./data/SRR7973881.bamstats)

### Step 2: Merge BAM files

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
The summary statistics of the final, merged dataset of mapped reads can be found here:
- [shad.merged](./data/shad.merged.bamstats)

### Step 3: Call SNPs
Variants were called with `bcftools v1.10.2 mpileup`, filtered, and output as a gzipped VCF.
```bash
# Call SNPs
threads=48 # CPU threads to use
bcftools \
        mpileup \
        -C 50 \
        -q 20 \
        -Q 20 \
        --threads ${threads} \
        -Ou \
        --ignore-RG \
        --per-sample-mF \
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
        -f Asap.fa \
        shad.merged.bam | \
        bcftools \
           call \
           --threads ${threads} \
           --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \
           -mv -Ou | \
        bcftools \
           filter \
           --threads ${threads} \
           -sFAIL -e'QUAL < 30 || INFO/MQ <= 30 || INFO/MQBZ < -9 || INFO/RPBZ < -5 || INFO/RPBZ > 5 || INFO/BQBZ < -5 || INFO/DP4[3]+INFO/DP4[4] <= 2' \
           -g10 \
           -G10 \
           -Ov | \
        bgzip -c > shad.vcf.gz
```

Next, to be compatible later with submission to the European Variation Archive, we had to rename the sample to match that in the BioSample database, Asap-001.  Also, to avoid any downstream issues with the mitochondrial contig, we removed all variants on the mitochondrial sequence [NC_014690.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_014690.1/).
```bash
# Rename the sample and remove mitochondrial contigs
   # The file "samples":
   # shad.merged.bam Asap-001
bcftools reheader \
	-s samples \
	shad.vcf.gz | \
	bcftools view \
		--targets ^NC_014690.1 \
		-O z \
		-o Asap001.noMito.vcf.gz -

# Index the VCF file
bcftools index --threads ${threads} Asap001.noMito.vcf.gz

# Stats
bcftools stats Asap001.noMito.vcf.gz > Asap001.noMito.stats

# Filter
bcftools view -f 'PASS' -O v Asap001.noMito.vcf.gz | bgzip -c > Asap001.filtered.noMito.vcf.gz
bcftools stats -f "PASS" Asap001.filtered.noMito.vcf.gz > Asap001.filtered.noMito.stats
```

The summary of SNP-calling statistics can be found here:
- Before Filtering: [Asap001.noMito.stats](./data/Asap001.noMito.stats)
- After Filtering: [Asap001.filtered.noMito.stats](./data/Asap001.filtered.noMito.stats)

### Step 4: Annotate SNPs
Variant information was done with `SnpEff v5.2e`.

```bash
# Build database
snpEff \
      build \
      -c snpEff.config \
      -gtf22 \
      -v -noCheckCds \
      -noCheckProtein \
      -nodownload \
      aloSap

# Annotate vcf
java -jar snpEff.jar ann \
	-c snpEff.config \
	aloSap \
	Asap001.filtered.noMito.vcf.gz | \
	bgzip -c > Asap001.filtered.ann.noMito.vcf.gz
```
The results from SnpEff can be found here:
- HTML format (open in web browser): [snpEff_summary.html](./data/snpEff_summary.html)
- PDF format: [snpEff_summary.pdf](./data/snpEff_summary.pdf)
- variant counts by gene (tab-delimited table): [snpEff_genes.txt](./data/snpEff_genes.txt)
