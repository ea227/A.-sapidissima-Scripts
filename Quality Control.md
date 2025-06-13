# **Quality Control**
### Step 1: Process paired-end reads

The paired-end reads were quality controlled and processed (including adapter removal) using [*fastp v0.20.0*](https://github.com/OpenGene/fastp) ([Chen et al. 2018](https://doi.org/10.1093/bioinformatics/bty560)). The details of various trimming and cleaning parameters are detailed below the code.

```bash
SRR=SRR7973879
tjreads=48

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
    --json=${out}.json \ 
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
  - ***--json=${SRR}.json*** :: output file name, JSON format
  - ***--html=${SRR}.html*** :: output file name, HTML format
  - ***--report_title="$SRR"*** :: output report tile
  - ***--thread=${threads}***  :: number of cpus to use

The output summary files from the read cleaning using _fastp_ can be found here (download and open the .html files in a web broser):
- html: [Pool B1: Early Spring](./data/B1.html) or as PDF [Pool B1: Early Spring](./data/B1.pdf)
- html: [Pool B2: Late Spring](./data/B2.html) or as PDF [Pool B2: Late Spring](./data/B2.pdf)
- html: [Pool C1: Early Autumn](./data/C1.html) or as PDF [Pool C1: Early Autumn](./data/C1.pdf)
- html: [Pool C2: Late Autumn](./data/C2.html) or as PDF [Pool C2: Late Autumn](./data/C2.pdf)



### Step 2: Process mate-pair reads
Mate pair libraries were processed with [NxTrim v0.4.3](https://github.com/sequencing/NxTrim). NxTrim removes the Nextera Mate Pair junction adapters and categorise reads according to the orientation implied by the adapter location. This allows the resulting reads to be formatted and oriented similar to paired-end reads for mapping.

```bash
# Process with NxTrim (one one example shown)
SRR=SRR7973880
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



fastp parameters: 

--in1/--in2: input forward and reverse read files, recognizes gzip

--stdout: write to standard out for piping

--adapter_fasta file: a file of known Illumina adapters to trim

--cut_front: enable a 5' sliding window trimmer

--cut_tail: enable a 3' sliding window trimmer

--cut_window_size: window size for the trimming

--cut_mean_quality: mean base score across the window required, or else trim the last base

--qualified_quality_phred: minimum base quality score to keep

--average_qual: remove read if the average quality across all bases is < specified value (20)

--unqualified_percent_limit: Percent of bases allowed to be less than q (--qualified_quality_phred) in a read

--n_base_limit: if read number of N bases is > specified value, then the read pair is discarded

--length_required: minimum read length to keep after trimming

--low_complexity_filter: filter sequences with a low complexity

--complexity_threshold: threshold for sequence complexity filter

--overrepresentation_analysis: look for overrepresented sequences, like adapters

--trim_poly_x: trim strings of homopolymers at the 3' end of reads

--poly_x_min_len: minimum length of homopolymer to trim

--json=${out}.json: output file name, JSON format

--html=${out}.html: output file name, HTML format

--report_title="$out": output report tile

--thread=${threads}: number of CPUs to use




