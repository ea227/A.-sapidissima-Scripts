
## **Data Accession** 

The short insert library and both mate pair libaries were accessed from NCBI with SRA toolkit v3.1.1:
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
fasterqdump {$SRR} --split files
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





