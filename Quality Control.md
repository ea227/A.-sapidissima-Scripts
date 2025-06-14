# **Quality Control**
### Step 1: Process paired-end reads

The paired-end reads were quality controlled and processed (including adapter removal) using [*fastp v0.20.0*](https://github.com/OpenGene/fastp) ([Chen et al. 2018](https://doi.org/10.1093/bioinformatics/bty560)). The details of various trimming and cleaning parameters are detailed below the code. A file of Illumina adapters to be targeted for removal was included and can be found here: [adapters.fa](./data/adapters.fa).

```bash
SRR=SRR7973879
threads=48

# Setup file of adapters to trim
adapters=adapters.fa

# Clean the reads using fastp
fastp \
    --in1=${SRR}_1.fastq.gz \ 
    --in2=${SRR}_2.fastq.gz \
    --adapter_fasta=${adapters} \ 
    --cut_front \ 
    --cut_tail \ 
    --cut_window_size=4 \ 
    --cut_mean_quality=20 \ 
    --qualified_quality_phred=20 \ 
    --average_qual=20 \ 
    --unqualified_percent_limit=30 \ 
    --n_base_limit=5 \ 
    --length_required=50 \ 
    --low_complexity_filter \ 
    --complexity_threshold=30 \ 
    --overrepresentation_analysis \ 
    --trim_poly_x \ 
    --poly_x_min_len=10 \ 
    --html=${out}.html \  
    --report_title="${SRR}" \ 
    --thread=${threads} \ 
    --out1 ${SRR}.finalclean_1.fq.gz \ 
    --out2 ${SRR}.finalclean_2.fq.gz
```
- _Parameters Explained:_
  - ***--in1/--in2*** :: input forward and reverse read files, recognizes gzip
  - ***--adapter_fasta file*** :: a file of known Illumina adapters to trim
  - ***--cut_front*** :: enable a 5' sliding window trimmer, like trimmomatic
  - ***--cut_tail*** :: enable a 3' sliding window trimmer, like trimmomatic
  - ***--cut_window_size=4*** :: window size for the trimming
  - ***--cut_mean_quality=20*** :: mean base score across the window required, or else trim the last base
  - ***--qualified_quality_phred=20*** :: minimum base quality score to keep
  - ***--average_qual=20*** :: remove read of the average quality across all bases is < 20
  - ***--unqualified_percent_limit=30*** :: Percent of bases allowed to be less than q in a read
  - ***--n_base_limit=5*** :: if one read's number of N bases is >5, then this read pair is discarded
  - ***--length_required=50*** :: minimum read length to keep after trimming
  - ***--low_complexity_filter*** :: filter sequences with a low complexity
  - ***--complexity_threshold=30*** :: threshold for sequence complexity filter
  - ***--overrepresentation_analysis*** :: look for overrepresented sequences, like adapters
  - ***--trim_poly_x*** :: trim strings of homopolymers at the 3' end of reads
  - ***--poly_x_min_len 10*** :: minimum length of homopolymer ot trim
  - ***--html=${SRR}.html*** :: output file name, HTML format
  - ***--report_title="$SRR"*** :: output report tile
  - ***--thread=${threads}***  :: number of cpus to use

The output summary files from the read cleaning using _fastp_ can be found here (download and open the .html files in a web broser):
- html: [SRR7973879.html](./data/SRR7973879.html) or as PDF [SRR7973879.pdf](./data/SRR7973879.pdf)
- html: [SRR7973880.html](./data/SRR7973880.html) or as PDF [SRR7973880.pdf](./data/SRR7973880.pdf)
- html: [SRR7973881.html](./data/SRR7973881.html) or as PDF [SRR7973881.pdf](./data/SRR7973881.pdf)

### Step 2: Process mate-pair reads
Mate pair libraries were processed with [NxTrim v0.4.3](https://github.com/sequencing/NxTrim) ([O'Connell et al. 2015](https://doi.org/10.1093/bioinformatics/btv057)). NxTrim removes the Nextera Mate Pair junction adapters and categorise reads according to the orientation implied by the adapter location. This allows the resulting reads to be formatted and oriented similar to paired-end reads for mapping. The output reads from _NxTrim_ were quality controlled and processed (including adapter removal) using [*fastp v0.20.0*](https://github.com/OpenGene/fastp) ([Chen et al. 2018](https://doi.org/10.1093/bioinformatics/bty560)) identical to the paired-end reads above.

```bash
# Process with NxTrim (one one example shown)
SRR=SRR7973880
 nxtrim \  #processing for mate-pair libraries
        -1 ${SRR}_1.fastq.gz \ 
        -2 ${SRR}_2.fastq.gz \ 
        -a \ 
        -l 50 \ 
        --stdout | \
       fastp \
              --stdin \
              --interleaved_in \
              --adapter_fasta=${adapters} \
              --cut_front \
              --cut_tail \
              --cut_window_size=4 \
              --cut_mean_quality=20 \
              --qualified_quality_phred=20 \
              --average_qual=20 \
              --unqualified_percent_limit=30 \
              --n_base_limit=5 \
              --length_required=50 \
              --low_complexity_filter \
              --complexity_threshold=30 \
              --overrepresentation_analysis \
              --trim_poly_x \
              --poly_x_min_len=10 \
              --html=${SRR}.html \
              --report_title="${SRR}" \
              --thread=${threads} \
              --out1 ${SRR}.finalclean_1.fq.gz \
              --out2 ${SRR}.finalclean_2.fq.gz
```

_NxTrim Parameters Explained:_
  - ***-1/-2*** :: input reads
  - ***--a*** :: aggressive adaptor search
  - ***-l*** :: minimum length of 50 bp after processing
  - ***--stdout*** :: pipe directly to standard out (directly to fastp)
