#!/usr/bin/python

import sys
import re
import subprocess
import os

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


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore.bed"):
    bed = re.split('[\t\n]',c)
    if float(bed[9]) >= 0:
        outfile.writelines(c)

outfile.close()



outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore.bed"):
    bed = re.split('[\t\n]',c)
    if float(bed[9]) >= 0:
        outfile.writelines(c)

outfile.close()


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore.bed"):
    bed = re.split('[\t\n]',c)
    if float(bed[9]) >= 0:
        outfile.writelines(c)

outfile.close()



outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore.bed"):
    bed = re.split('[\t\n]',c)
    if float(bed[9]) >= 0:
        outfile.writelines(c)

outfile.close()

