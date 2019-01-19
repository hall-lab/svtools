import argparse, sys, StringIO
import pandas as pd
import numpy as np
import scipy.spatial.distance as ssd
import pysam
sys.path.insert(1,'/gscmnt/gc2802/halllab/abelhj/svtools')
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from collections import namedtuple
import svtools.utils as su



def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--vcf', metavar='<VCF>', dest='manta_vcf', help="manta input vcf")
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true')
    parser.add_argument('-S', '--slop', dest='slop',  default=0,  required=False, help='padding to either side')

def command_parser():
    parser = argparse.ArgumentParser(description="cross-cohort cnv caller")
    add_arguments_to_parser(parser)
    return parser

def convert_variant(v):
    set_read_counts(v)
    if v.get_info('SVTYPE')=='DEL':
        convert_del(v)
    elif v.get_info('SVTYPE')=='DUP':
        convert_dup(v)
    elif v.get_info('SVTYPE')=='INV':
        convert_inv(v)
    elif v.get_info('SVTYPE')=='BND':
        convert_bnd(v)

def split_ci(ci):
     return[int(ci.split(',')[0]),  int(ci.split(',')[1])]

def uniform_pr(length):
    pr=np.ones(length, dtype='float64')/length
    pr1=','.join( map(str, pr))
    return pr1

def set_read_counts(var):

    sample=var.sample_list[0]
    gt=var.genotype(sample)
    pe=0
    sr=0
    if 'PR' in var.format_dict:
        pe=int(gt.get_format('PR').split(',')[1])
    if 'SR' in var.format_dict:
        sr=int(gt.get_format('SR').split(',')[1])
    var.info['PE']=pe
    var.info['SR']=sr
    var.info['SU']=pe+sr

def set_cis_prs(v):
    imprec=False
    cipos='0,0'
    ciend='0,0'
    prpos=1.0
    prend=1.0
    if 'CIPOS' in v.info:
        cipos=v.get_info('CIPOS')
        [start, stop]=split_ci(cipos)
        prpos=uniform_pr(stop-start+1)
        imprec=True
    if 'CIEND' in v.info:
        ciend=v.get_info('CIEND')
        [start, stop]=split_ci(ciend)
        prend=uniform_pr(stop-start+1)
        imprec=True
    v.info['CIPOS']=cipos
    v.info['CIEND']=ciend
    v.info['PRPOS']=prpos
    v.info['PREND']=prend
    v.set_info('IMPRECISE', imprec)
    
def convert_del(var):
    var.alt='<DEL>'
    var.info['STRANDS']="+-:1"
    var.ref='N'

def convert_dup(var):
    var.alt='<DUP>'
    var.info['STRANDS']="-+:1"
    var.ref='N'

def convert_inv(var):
    var.ref='N'
    var.alt='<INV>'
    if 'INV3' in var.info:
        var.info['STRANDS']="++:1"
    else:
        var.info['STRANDS']="--:1"

def convert_bnd(var):
    var.ref='N'
    alt=var.alt
    ff=alt.find("[")
    newalt=""
    strands=""
    if ff==0:
        strands="--:6"
        ff1=alt.find("[", 1)
        newalt=alt[0:(ff1+1)]+'N'
    elif ff>0:
        strands="+-:6"
        newalt='N'+alt[ff::]
    else:
        ff=alt.find("]")
        if ff==0:
            strands="-+:6"
            ff1=alt.find("]", 1)
            newalt=alt[0:(ff1+1)]+'N'
        else:
            strands="++:6"
            newalt='N'+alt[ff::]
    print(var.alt+"\t"+newalt+"\t"+strands)
    var.alt=newalt
    var.info['STRANDS']=strands
        

def run_from_args(args):

  vcf = Vcf()
  vcf_out=sys.stdout
  in_header = True
  header_lines = list()
  with su.InputStream(args.manta_vcf) as input_stream:
    for line in input_stream:
      if in_header:
        header_lines.append(line)
        if line[0:6] == '#CHROM':
          in_header=False
          vcf.add_header(header_lines)
          vcf.add_info('PRPOS', '1', 'String', 'Breakpoint probability dist')
          vcf.add_info('PREND', '1', 'String', 'Breakpoint probability dist')
          vcf.add_info('STRANDS', '.', 'String', 'Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--')
          vcf.add_info('SU', '.', 'Integer', 'Number of pieces of evidence supporting the variant across all samples')
          vcf.add_info('PE', '.', 'Integer', 'Number of paired-end reads supporting the variant across all samples')
          vcf.add_info('SR', '.', 'Integer', 'Number of split reads supporting the variant across all samples')
          vcf_out.write(vcf.get_header()+'\n')
      else:
        v = Variant(line.rstrip().split('\t'), vcf)
        convert_variant(v)
        vcf_out.write(v.get_var_string()+"\n")
          

parser=command_parser()
args=parser.parse_args()
print(str(args))
run_from_args(args)
