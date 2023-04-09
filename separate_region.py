#!/usr/bin/python

import sys
import re
import subprocess
import os
import numpy as np
from numba import jit
import pandas as pd
from tqdm import tqdm
from time import sleep


#-----------------
#Usage or get argv
#-----------------
argvs = sys.argv
argc = len(argvs)
if (argc == 2):
    chr = argvs[1]
else:
    errmsg = "Usage:  " + argvs[0] + "  ID_for_a_individual"
    print ('')
    print (errmsg)
    print ('')
    sys.exit()





print (chr)
#------------------------
#Exon region slop 50 nt
#------------------------
#Make chromesome size file
#cut -f 1,2 /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa.fai > /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.chrom.sizes

cmd = "bedtools slop -i /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + ".gtf -g /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.chrom.sizes -b 50 > /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + "_slop50.gtf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)
#------------------
#Acceptor site (A)
#------------------
print ("accA")
cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + ".header /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01.vcf > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h.vcf -b /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + "_slop50.gtf -v |sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_deepintron.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h.vcf -b /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + "_slop50.gtf -wa |sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)


#------------------
#Acceptor site (G)
#------------------
print ("accG")
cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + ".header /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01.vcf > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h.vcf -b /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + "_slop50.gtf -v |sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_deepintron.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h.vcf -b /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + "_slop50.gtf -wa |sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

#------------------
#Donor site (G)
#------------------
print ("donoG")

cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + ".header /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01.vcf > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h.vcf -b /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + "_slop50.gtf -v |sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_deepintron.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h.vcf -b /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + "_slop50.gtf -wa |sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)



#------------------
#Donor site (T)
#------------------
print ("donoT")

cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + ".header /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01.vcf > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h.vcf -b /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + "_slop50.gtf -v |sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_deepintron.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h.vcf -b /mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + "_slop50.gtf -wa |sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)
