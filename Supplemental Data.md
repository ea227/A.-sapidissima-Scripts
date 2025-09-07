## **Variant Density Across the Genome**

Variant information was accessed for 20000 bp windows using vcftools v0.1.17, and output as a .tsv file: 
```
vcftools --gzvcf Asap001.filtered.ann.noMito.vcf.gz --SNPDensity 20000 --out Asap001
```
The following plot shows variant density visualized in windows of 20000 bp, separated by chromosome number: 
[variantdensityplot.pdf](https://github.com/user-attachments/files/22197407/variantdensityplot.pdf)


#old
[variantdensityplot.pdf](https://github.com/user-attachments/files/19489152/variantdensityplot.pdf)

## **fastp Reports**

Short insert library:
[srr79.pdf](https://github.com/user-attachments/files/19489292/srr79.pdf)

Mate pair library (5-7kb insert):
[srr80.pdf](https://github.com/user-attachments/files/19489294/srr80.pdf)

Mate pair library (10-12kb insert):
[srr81.pdf](https://github.com/user-attachments/files/19489296/srr81.pdf)




