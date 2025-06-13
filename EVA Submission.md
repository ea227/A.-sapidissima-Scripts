# Submitting the variant data to the EVA
The code below is for submitting the variant data (i.e., vcf file) to the European Variation Archive ([EVA](https://www.ebi.ac.uk/eva/))

## Step 1: Preparing VCF for EVA
I first had to change the sample name to match exactly between the VCF file and the BioSample info. I used bcftools reheader to do this using the following command:
```bash
# Rename the samples in the VCF
bcftools reheader \
   -s samples \
   -o Asap001.filtered.ann.vcf.gz \
   shad.filtered.ann.vcf.gz

# The file "samples":
shad.merged.bam Asap-001
```

2. I had to remove the mitochondrial contig:
```
# Remove mitochondrial contig NC_014690.1
bcftools view \
   --targets ^NC_014690.1 \
   -O z \
   -o Asap001.filtered.ann.noMito.vcf.gz \
   Asap001.filtered.ann.vcf.gz
```
3. I then had to remove the following line (was offending to the eva-submit software):
```
# Unzip the vcf file
gunzip Asap001.filtered.ann.noMito.vcf.gz

# Used vim to delete the following line and saved
##ALT=<ID=*,Description="Represents allele(s) other than observed.">

# Recompressed the vcf
bgzip Asap001.filtered.ann.noMito.vcf
```
4. Installed the eva submission software using conda
```
# Install eva-sub-cli
conda create \
   -n eva \
   -c conda-forge \
   -c bioconda \
   eva-sub-cli
conda activate eva

# Validate via eva-sub-cli
eva-sub-cli.py \
   --metadata_xlsx \
   Asap-001_EVA_Submission.V2.0.1.xlsx \
   --submission_dir . \
   --tasks VALIDATE \
   --debug

# Submit via eva-sub-cli
eva-sub-cli.py \
   --metadata_xlsx Asap-001_EVA_Submission.V2.0.1.xlsx \
   --submission_dir . \
   --tasks SUBMIT \
   --debug \
   --username Webin-70646 \
   --password '##########'
```

