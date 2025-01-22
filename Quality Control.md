Cleaned sequence data for each library was wuality controlled and preprocessed using fastp v0.20.0:
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
