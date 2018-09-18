# Tiered merging for large cohorts
For very large SV callsets, we recommend a tiered approach to merging the individual Lumpy VCFs. This is useful to keep compute requirements modest and also to help smooth out batch effects between cohorts. For our own large cohorts, we've adopted a merging strategy whereby we merge groups of up to 1000 samples per cohort (larger cohorts will have multiple batches of 1000 or less) and then sort and merge the subsequent merged files again.

## Initial Per-batch sorting and merging
This is easily done by providing a file input containing the path to each individual Lumpy VCF for each group of 1000 or less samples and then subsequently merging.

```
svtools lsort -f batch_of_lumpy_vcfs.txt 
  | svtools lmerge -i /dev/stdin -f 20 
  | bgzip -c > batch.merged.vcf.gz
```

## Final sorting and merging
After this step you will have one output file per batch. However, these files will _not_ contain genotypes so you'll need to specify additional options to ensure that they are properly combined.

```
svtools lsort -r -f file_of_merged_batches --batch-size 1
  | svtools lmerge -i /dev/stdin -f 20 -w carrier_wt
  | bgzip -c > final_output.merged.vcf.gz
```
