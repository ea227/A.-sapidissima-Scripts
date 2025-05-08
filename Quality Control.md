## **Quality Control**

Cleaned sequence data for each library was quality controlled and preprocessed using fastp v0.20.0. Output html files are also included as PDFs in the Supplemental Data section. The following parameters were used:
```
fastp \ #Quality control and preprocessing of raw sequencing data  
    --in1=${r1} \ 
    --in2=${r2} \ 
    --stdout \ 
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
    --report_title="${out}" \ 
    --thread=${threads} \ 
    --out1 cleaned/${SRR}.finalclean_1.fq.gz \ 
    --out2 cleaned/${SRR}.finalclean_2.fq.gz
```
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




