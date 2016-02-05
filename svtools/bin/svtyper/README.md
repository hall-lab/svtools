svtyper
=======

Bayesian genotyper for structural variants

### Example workflow

#### Data
```
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.bam
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.bam.bai
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.splitters.bam
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.splitters.bam.bai
wget http://colbychiang.com/hall/cshl_sv_2014/data/NA12878.20.vcf.gz
```

#### Genotype with SVTyper
```
zcat NA12878.20.vcf.gz \
    | ./svtyper \
        -B NA12878.20.bam \
        -S NA12878.20.splitters.bam \
        > NA12878.20.gt.vcf
```
#### Warning
2015-10-05

As of commit [2c2ef7f91698a6d2929430f0865402ad421a8e3d](https://github.com/hall-lab/svtyper/commit/2c2ef7f91698a6d2929430f0865402ad421a8e3d), SVTyper assumes that BAMs were aligned with BWA MEM **without** the "-M" flag. If you used the "-M" flag in your alignment, then you should also use the "-M" flag when running SVTyper.
