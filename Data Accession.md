
## **Data Accession** 

The short insert library and both mate pair libaries were accessed from NCBI with SRA toolkit v3.1.1. {SRR} was substituted for the SRA accession numbers for each library:
```
prefetch ${SRR}
```

The genome was accessed and indexed using bwa v0.7.18:
```
bwa index ${SRR}
```
