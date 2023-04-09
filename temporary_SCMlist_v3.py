#!/usr/bin/python

import sys
import re
import subprocess
import os
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

#----------------------
#Acceptor site (A)
#----------------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_snv.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0.bed"):
    bed =  re.split('[\t\n]',c)
    if "+" == str(bed[5]):
        snv = int(bed[1]) + 19
        outfile.writelines(c[:-1] + "\t" + str(snv) + "\n")
    if "-" == str(bed[5]):
        snv = int(bed[1]) + 5
        outfile.writelines(c[:-1] + "\t" + str(snv) + "\n")
outfile.close()


#----------------------
#Acceptor site (G)
#----------------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_snv.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0.bed"):
    bed =  re.split('[\t\n]',c)
    if "+" == str(bed[5]):
        snv = int(bed[1]) + 20
        outfile.writelines(c[:-1] + "\t" + str(snv) + "\n") 
    if "-" == str(bed[5]):
        snv = int(bed[1]) + 4
        outfile.writelines(c[:-1] + "\t" + str(snv) + "\n")
outfile.close()
    

#----------------------
#Donor site (G)
#----------------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_snv.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0.bed"):
    bed =  re.split('[\t\n]',c)
    if "+" == str(bed[5]):
        snv = int(bed[1]) + 4
        outfile.writelines(c[:-1] + "\t" + str(snv) + "\n")
    if "-" == str(bed[5]):
        snv = int(bed[1]) + 6
        outfile.writelines(c[:-1] + "\t" + str(snv) + "\n")
outfile.close()
    
#----------------------
#Donor site (T)
#----------------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_snv.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0.bed"):
    bed =  re.split('[\t\n]',c)
    if "+" == str(bed[5]):
        snv = int(bed[1]) + 5
        outfile.writelines(c[:-1] + "\t" + str(snv) + "\n")
    if "-" == str(bed[5]):
        snv = int(bed[1]) + 5
        outfile.writelines(c[:-1] + "\t" + str(snv) + "\n")
outfile.close()
    
