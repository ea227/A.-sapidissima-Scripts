**Introduction**

The bioinformatic analysis of the raw sequencing reads described here, including links and citations to programs, and software versions, and
the scripts used. All analyses were run on the coombs computer cluster at the University of Central Florida (UCF).
Coombs is a Linux server currently running Ubuntu 20.04.6 LTS.

All sequence data referenced herein are attributed to an award from the The Sequencing and Genomic Technologies Shared Resource at the Duke 
Center for Genomic and Computational Biology (Duke University). Sequencing was done by (R Fitak, etc. edit this). Samples and assistance were provided 
by Stephen Jackson and Kevin Dockendorf at the Edenton National Fish Hatchery. (bob edit?). 

The reference genome was sequenced and assembled at the Vertebrate Genomes Lab at Rockefeller University for the Vertebrate Genomes Project (VGP) and coordinated by Olivier Fedrigo and Erich D. Jarvis.
Reid Hyle collected a sample for Nina Overgaard Therkildsen from the St. Johns River, Florida, USA for genome assembly. 

**Data Accession** 

The short insert library and both mate pair libaries were accessed from NCBI with SRA toolkit v3.1.1:
```
prefetch ${SRR}
```

The genome was accessed and indexed using bwa v0.7.18:
```
bwa index ${SRR}
```

**Data Cleaning**

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





