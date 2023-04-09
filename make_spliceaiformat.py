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

#------------------
#Acceptor site (A)
#------------------
print ("accA")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_snv.vcf","w")

f = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.vcf")
pbar = tqdm(total=len(open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.vcf").readlines()))
for c in tqdm(f):
    pbar.update(1)
    vcf =  re.split('[\t\n]',c)
    #print (vcf[-2])
    for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_snv.bed"):
        bed = re.split('[\t\n]',i)
        #print (str(vcf[0]),str(bed[0]),str(vcf[1]),str(bed[10]),str(vcf[-2]), str(bed[5]),str(vcf[3]),str(bed[3]),str(vcf[4]),str(bed[4]))
        if str(vcf[0]) == str(bed[0]) and str(vcf[1]) == str(bed[10]) and str(vcf[-2]) == str(bed[5]) and str(vcf[3]) == str(bed[3]) and str(vcf[4]) == str(bed[4]):
            outfile.writelines(c)
            break
outfile.close()
pbar.close()



#------------------
#Acceptor site (G) and Donor site (G)
#------------------
print ("donoaccG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_snv.vcf","w")
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_snv.vcf","w")
f = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoaccG.vcf")
pbar = tqdm(total=len(open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoaccG.vcf").readlines()))
for c in tqdm(f):
    pbar.update(1)
    vcf =  re.split('[\t\n]',c)
    #print (vcf[-2])
    for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_snv.bed"):
        bed = re.split('[\t\n]',i)
        #print (str(vcf[0]),str(bed[0]),str(vcf[1]),str(bed[10]),str(vcf[-2]), str(bed[5]),str(vcf[3]),str(bed[3]),str(vcf[4]),str(bed[4]))
        if str(vcf[0]) == str(bed[0]) and str(vcf[1]) == str(bed[10]) and str(vcf[-2]) == str(bed[5]) and str(vcf[3]) == str(bed[3]) and str(vcf[4]) == str(bed[4]):
            outfile.writelines(c)
            break
    for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_snv.bed"):
        bed = re.split('[\t\n]',i)
        #print (str(vcf[0]),str(bed[0]),str(vcf[1]),str(bed[10]),str(vcf[-2]), str(bed[5]),str(vcf[3]),str(bed[3]),str(vcf[4]),str(bed[4]))
        if str(vcf[0]) == str(bed[0]) and str(vcf[1]) == str(bed[10]) and str(vcf[-2]) == str(bed[5]) and str(vcf[3]) == str(bed[3]) and str(vcf[4]) == str(bed[4]):
            outfile2.writelines(c)
            break

    
outfile.close()
outfile2.close()
pbar.close()


#------------------
#Donor site (T)
#------------------
print ("donoT")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_snv.vcf","w")

f = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.vcf")
pbar = tqdm(total=len(open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.vcf").readlines()))
for c in tqdm(f):
    pbar.update(1)
    vcf =  re.split('[\t\n]',c)
    #print (vcf[-2])
    for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_snv.bed"):
        bed = re.split('[\t\n]',i)
        #print (str(vcf[0]),str(bed[0]),str(vcf[1]),str(bed[10]),str(vcf[-2]), str(bed[5]),str(vcf[3]),str(bed[3]),str(vcf[4]),str(bed[4]))
        if str(vcf[0]) == str(bed[0]) and str(vcf[1]) == str(bed[10]) and str(vcf[-2]) == str(bed[5]) and str(vcf[3]) == str(bed[3]) and str(vcf[4]) == str(bed[4]):
            outfile.writelines(c)
            break
outfile.close()
pbar.close()


