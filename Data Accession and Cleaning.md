** Introduction **

The bioinformatic analysis of the raw sequencing reads described here, including links and citations to programs, and software versions, and
the scripts used. All analyses were run on the coombs computer cluster at the University of Central Florida (UCF).
Coombs is a Linux server currently running Ubuntu 20.04.6 LTS.

All sequence data referenced herein are attributed to an award from the The Sequencing and Genomic Technologies Shared Resource at the Duke 
Center for Genomic and Computational Biology (Duke University). Sequencing was done by (R Fitak, etc. edit this). Samples and assistance were provided 
by Stephen Jackson and Kevin Dockendorf at the Edenton National Fish Hatchery. (bob edit?). 

The reference genome was sequenced and assembled at the Vertebrate Genomes Lab at Rockefeller University for the Vertebrate Genomes Project (VGP) and coordinated by Olivier Fedrigo and Erich D. Jarvis.
Reid Hyle collected a sample for Nina Overgaard Therkildsen from the St. Johns River, Florida, USA for genome assembly. 

**Data Accession** 
The short insert library and both mate pair libaries were accessed from NCBI by:
prefetch ${SRR}

The genome was accessed and indexed by:
bwa index ${SRR}




