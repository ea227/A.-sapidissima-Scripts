# **Data Accession and Preparation** 

### Step 1: Obtain raw sequencing reads
The paired-end (short insert) library and both mate-pair libaries were accessed from NCBI with SRA toolkit v3.1.1.
For the paired-end reads, `SRR=SRR7973879`, and for the mate-pair reads `SRR=SRR7973880` and `SRR=SRR7973881`.

```bash
# Download reads (paired-end reads example shown)
threads=48
SRR=SRR7973879
fasterq-dump \
        --progress \
        --threads ${threads} \
        --split-files \
        ${SRR}
```
_Parameters explained:_
- _--progress_ :: show progess bar during download
- _--threads_ :: the number of cpu threads to use, increases parallelization
- --split_files :: splits the fastq read files into separate forward (${SRR}_1.fastq.gz) and reverse (${SRR}_2.fastq.gz) read files.

Alternatively, there were occassional download issues using `fasterq-dump`. If so, then read files were first downoaded using `prefetch` in the `.sra` format and then converted to fastq format like this:
```bash
prefetch --max-size 200G -o ${SRR}.sra ${SRR}
threads=48
fasterq-dump \
        --progress \
        --threads ${threads} \
        --split-files \
        ./${SRR}.sra
```

### Step 2: Download and index the reference genome
The American shad reference genome (NCBI accession [GCF_018492685.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018492685.1/)) was downloaded and indexed using `bwa v0.7.18`:

```bash
# Download and index reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/492/685/GCF_018492685.1_fAloSap1.pri/GCF_018492685.1_fAloSap1.pri_genomic.fna.gz
gunzip GCF_018492685.1_fAloSap1.pri_genomic.fna.gz
mv GCF_018492685.1_fAloSap1.pri_genomic.fna Asap.fa
bwa index Asap.fa
```
