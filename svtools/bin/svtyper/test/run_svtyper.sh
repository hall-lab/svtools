# zcat NA12878.20.vcf.gz \
#     | ../svtyper \
#     -n 10000 \
#     -B /gscmnt/gc2719/halllab/users/cchiang/data/toy_data/NA12878.gtex.realign.chrom20/NA12878.gtex.realign.20.bam \
#     -S /gscmnt/gc2719/halllab/users/cchiang/data/toy_data/NA12878.gtex.realign.chrom20/NA12878.gtex.realign.20.splitters.bam


SVTYPER=../svtyper
# SVTYPER=/gscmnt/gc2719/halllab/users/cchiang/src/svtyper/svtyper

cat NA12893.het.vcf \
    | $SVTYPER \
    -B /gscmnt/gc2802/halllab/ceph1463_realign_021815/NA12893/NA12893.bam \
    -S /gscmnt/gc2802/halllab/ceph1463_realign_021815/NA12893/NA12893.bam \
    > NA12893.het.gt.vcf

echo "Heterozygous computed genotypes:"
cat NA12893.het.gt.vcf | vawk '{ print S$NA12893$GT }' | sort | uniq -c

for TYPE in DEL DUP INV
do
    cat NA12893.het.gt.vcf | vawk -v TYPE=$TYPE "{ if (I\$SVTYPE==TYPE && I\$SVLEN>=-1000000 && I\$SVLEN<=1000000) print S\$NA12893\$GT }" | sort | uniq -c | awk -v TYPE=$TYPE 'BEGIN { MATCH=0; INF=0 } { if ($2=="0/1") MATCH+=$1 ; INF+=$1 } END { print TYPE,MATCH,INF,MATCH/INF }' OFS="\t"
done

echo -e "\nHomozygous computed genotypes:"
cat NA12893.hom.vcf \
    | $SVTYPER \
    -B /gscmnt/gc2802/halllab/ceph1463_realign_021815/NA12893/NA12893.bam \
    -S /gscmnt/gc2802/halllab/ceph1463_realign_021815/NA12893/NA12893.bam \
    > NA12893.hom.gt.vcf

cat NA12893.hom.gt.vcf | vawk '{ print S$NA12893$GT }' | sort | uniq -c

for TYPE in DEL DUP INV
do
    cat NA12893.hom.gt.vcf | vawk -v TYPE=$TYPE "{ if (I\$SVTYPE==TYPE && I\$SVLEN>=-1000000 && I\$SVLEN<=1000000) print S\$NA12893\$GT }" | sort | uniq -c | awk -v TYPE=$TYPE 'BEGIN { MATCH=0; INF=0 } { if ($2=="1/1") MATCH+=$1 ; INF+=$1 } END { print TYPE,MATCH,INF,MATCH/INF }' OFS="\t"
done
