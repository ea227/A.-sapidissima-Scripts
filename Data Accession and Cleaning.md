## **Data Accession** 
The short insert library and both mate-pair libaries were accessed from NCBI with SRA toolkit v3.1.1.
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
_--progress_ :: show progess bar during download
_--threads_ :: the number of cpu threads to use, increases parallelization
--split_files :: splits the fastq read files into separate forward (${SRR}_1.fastq.gz) and reverse (${SRR}_2.fastq.gz) read files.

Alternatively, there were occassional download issues using `fasterq-dump`. If so, then read files were first downoaded using `prefetch` in the `.sra` format and then converted to fastq format like this:
```bash
prefetch --max-size 200G -o ${SRR}.sra ${SRR}
fasterq-dump \
        --progress \
        --threads ${threads} \
        --split-files \
        ./${SRR}.sra
```








{SRR} was substituted for the SRA accession numbers for each library:
```
prefetch ${SRR}
```

The genome was accessed and indexed using bwa v0.7.18:
```
bwa index ${SRR}
```
## **Data Cleaning**

Short read libraries were converted to fastq format using fasterqdump from SRA tools:
```
fasterqdump ${SRR} --split files
```
Mate pair libraries were converted to fastq format using fastq-dump from SRA tools and then processed with NxTrim v0.4.3:
```
fastq-dump --gzip --split-files ${SRR} #Convert mate pair libraries to fastq 
 nxtrim \  #processing for mate-pair libraries
        -1 ${SRR}_1.fastq.gz \ 
        -2 ${SRR}_2.fastq.gz \ 
        -a \ 
        -l 50 \ 
        --stdout
```
NxTrim parameters:

--a: aggressive adaptor search

--l: minimum length 

--stdout: pipe directly to aligner 


