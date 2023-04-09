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



#------------
#all cat
#------------

print ("all")
cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + ".header /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "spliceai -I /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.vcf -O /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf -R /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa -A /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/GRCh38_gtf/v29/gencode_v29_spliceai.txt"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

#---------------
#acceptor A
#----------------

print ("accA")


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf"):
	vcf = re.split('[\t\n]',i)
	for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf"):
		if "#" in c:
			pass
		else:
			spliceai = re.split('[\t\n]',c)
			if str(vcf[0]) == str(spliceai[0]) and str(vcf[1]) == str(spliceai[1]) and str(vcf[2]) == str(spliceai[2]) and str(vcf[3]) == str(spliceai[3]) and str(vcf[4]) == str(spliceai[4]):
				outfile.writelines(c)
				break
outfile.close()

#---------------
#acceptor G
#----------------

print ("accG")

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf"):
        vcf = re.split('[\t\n]',i)
        for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf"):
                if "#" in c:
                        pass
                else:
                        spliceai = re.split('[\t\n]',c)
                        if str(vcf[0]) == str(spliceai[0]) and str(vcf[1]) == str(spliceai[1]) and str(vcf[2]) == str(spliceai[2]) and str(vcf[3]) == str(spliceai[3]) and str(vcf[4]) == str(spliceai[4]):
                                outfile.writelines(c)
                                break
outfile.close()

#---------------
#donor G
#----------------

print ("donoG")

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf"):
        vcf = re.split('[\t\n]',i)
        for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf"):
                if "#" in c:
                        pass
                else:
                        spliceai = re.split('[\t\n]',c)
                        if str(vcf[0]) == str(spliceai[0]) and str(vcf[1]) == str(spliceai[1]) and str(vcf[2]) == str(spliceai[2]) and str(vcf[3]) == str(spliceai[3]) and str(vcf[4]) == str(spliceai[4]):
                                outfile.writelines(c)
                                break
outfile.close()

#---------------
#donor T
#----------------

print ("donoT")

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50.vcf"):
        vcf = re.split('[\t\n]',i)
        for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf"):
                if "#" in c:
                        pass
                else:
                        spliceai = re.split('[\t\n]',c)
                        if str(vcf[0]) == str(spliceai[0]) and str(vcf[1]) == str(spliceai[1]) and str(vcf[2]) == str(spliceai[2]) and str(vcf[3]) == str(spliceai[3]) and str(vcf[4]) == str(spliceai[4]):
                                outfile.writelines(c)
                                break
outfile.close()
