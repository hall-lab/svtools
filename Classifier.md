# Variant reclassification with `svtools classify`
Variant reclassification identifies low quality deletions and duplications and changes their type to `BND`. The method by which this is down varies depending on sample size.
The classifier can be run in several modes depending on the sample size. For cohorts with >30 samples we recommend using the 'large_sample' mode. For smaller cohorts, we recommend the 'naive_bayes' mode. An experimental 'hybrid' mode that combines the two modes is also available, but not yet recommended.

## Running the classifier in 'large_sample' mode

In 'large_sample' mode, copynumber estimates are regressed against allele balance (for high frequency variants) or compared between individuals that are homozygous reference and non-reference (for rare variants). This is the preferred method if you have a large number of samples in your cohort.

```
zcat merged.sv.new_pruned.vcf.gz \
|  svtools classify \
 -g ceph.sex.txt \
 -a repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.b37.sorted.bed.gz \
 -m large_sample \
| bgzip -c > output.ls.vcf.gz
```

## Creating a set of high-quality, simple deletions and duplications for classifier training
`svtools classify` in the naive bayes mode requires a training set of high-quality deletions and duplications. A set of these can be identified as described below using the `mean_cn.pl` located in the scripts directory of the `svtools` repository. A VCF from a large cohort (>30 samples) should be used.

### Find the mean per-site copy number and overall percentiles for deletions and duplications
This step generates two files: one containing the mean copynumber of heterozygous and homozygous reference samples at all deletion and duplication sites, and one containing the 10th and 90th percentiles of these quantities over all sites.

```
zcat merged.sv.pruned.vcf.gz \
| vawk --header '{if((I$SVTYPE==“DEL” || I$SVTYPE==“DUP” || I$SVTYPE==“MEI”) && I$AF>0.01 && $1!="X" && $1!="Y") print $0}' \
| perl mean_cn.pl 1>per_site.means.txt 2>overall_percentiles.txt
```

### Extract high-quality duplication and deletion sites
This step extracts duplication and deletion sites where the the mean heterozygous copy number and mean homozygous reference copy number fall within the 10th and 90th percentiles of all samples. These sites are considered high-quality and can be intersected with a callset to be used as training data for `svtools classify`.
```
cat per_site.means.txt  \
| cut -f -8 \
| zjoin -a stdin -b <(cat overall_percentiles.txt | cut -f -8 ) -1 2 -2 1 \
| awk '{if($5>$11 && $5<$12 && $8>$15 && $8<$16) print $0}' \
| cut -f 1 \
| zjoin -a <(zcat sv.vcf.gz) -b stdin -1 3 -2 1 \
| svtools vcftobedpe \
| bgzip -c > training_vars.bedpe.gz
```
