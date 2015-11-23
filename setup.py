from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES

setup(
    name='svtools',
    version='0.1.1',
    author='Ira Hall lab',
    author_email='abadve@genome.wustl.edu',
    data_files=[('data', ['data/repeatMasker.recent.lt200millidiv.b37.sorted.bed.gz']),('data', ['data/repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.b37.sorted.bed.gz'])],
    scripts=["bin/vcftobedpe","bin/varlookup","bin/svtools","bin/vcfsort",\
             "bin/bedpesort","bin/prune","bin/copynumber","bin/bedpetovcf",\
             "bin/bedpetobed12","bin/afreq","bin/classify","bin/genotype",\
             "bin/lsort","bin/lmerge","bin/vcfpaste",\
             'scripts/l_bp.py',"scripts/vcf_group_multiline.py",\
             "bin/cnvnator-multi","scripts/create_coordinates.py","scripts/create.master.stats.sh","scripts/sv_counts.sh","scripts/get_stats_final.sh"],
    url='https://github.com/hall-lab/svtools',
    license='LICENSE.txt',
    description='Tools for processing and analyzing structural variants',
    long_description=open('README.txt').read()
)

#
# setup(
#     ...
#     dependency_links=['http://github.com/user/repo/tarball/master#egg=package-1.0']
#     ...
# )