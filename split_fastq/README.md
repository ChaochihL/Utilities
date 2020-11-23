# Alternative split fastq method

Script that splits large (>100GB gzipped) fastq files when Linux commands like `gzip -dc ${R1_FILE} | split --numeric-suffixes=1 -l ${NUM_LINES} - ${SAMPLE_PREFIX}_part --suffix-length=1 --additional-suffix=_R1.fastq` are producing some truncated sequences in the output when those sequences in the original gzipped fastq file are not truncated.

The script is called `split_fastq.py`. This subdirectory also contains some toy fastq files and a FastQC report for the toy fastq files.

Here is an example of how to run the script:

```bash
fastq_file="~/GitHub/Utilities/split_fastq/toy.fastq.gz"
num_parts="3" # number of split part files we want
total_num_seqs="10" # Obtained from FastQC report toy_fastqc.html
out_dir="~/Downloads"

./split_fastq.py ${fastq_file} \
    ${num_parts} \
    ${total_num_seqs} \
    ${out_dir}
```
