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

print ("accA")
#------------------
#Acceptor site (A)
#------------------
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf","w")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over05.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf"):
    if "#" in i:
        pass
    else:
        vcf = re.split('[\t\n;]',i)
        x = str(vcf[7])
        x = str(vcf[-2]) 
        for c in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_transcript_ID.txt"):
            IDs = re.split('[\t\n]',c)
            if str(IDs[0]) in x: 
                transcript = str(vcf[4]) + "|" + str(IDs[0])
                scores = re.split('[\t\n,]',x)
                for score in scores:
                    if str(transcript) in score:
                        SpliceAI = re.split('[\t\n|]',score)
                        value = float(SpliceAI[2])
                        if float(value) >= 0.50:
                            value = "{0:.2f}".format(float(value))
                            outfile.writelines(i[:-1] + "\t" + str(value) +  "\t" + str(IDs[0]) + "\n")
                        if float(value) >= 0.80:
                            value = "{0:.2f}".format(float(value))
                            outfile2.writelines(i[:-1] + "\t" + str(value) +  "\t" + str(IDs[0]) + "\n")

outfile.close()
outfile2.close()


#------------------
#Acceptor site (G)
#------------------
print ("accG")
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf","w")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over05.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf"):
    if "#" in i:
        pass
    else:
        vcf = re.split('[\t\n;]',i)
        x = str(vcf[7])
        x = str(vcf[-2]) 
        for c in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_transcript_ID.txt"):
            IDs = re.split('[\t\n]',c)
            if str(IDs[0]) in x: 
                transcript = str(vcf[4]) + "|" + str(IDs[0])
                scores = re.split('[\t\n,]',x)
                for score in scores:
                    if str(transcript) in score:
                        SpliceAI = re.split('[\t\n|]',score)
                        value = float(SpliceAI[2])
                        if float(value) >= 0.50:
                            value = "{0:.2f}".format(float(value))
                            outfile.writelines(i[:-1] + "\t" + str(value) +  "\t" + str(IDs[0]) + "\n")
                        if float(value) >= 0.80:
                            value = "{0:.2f}".format(float(value))
                            outfile2.writelines(i[:-1] + "\t" + str(value) +  "\t" + str(IDs[0]) + "\n")

outfile.close()
outfile2.close()

#------------------
#Donor site (G)
#------------------
print ("donoG")
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf","w")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over05.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf"):
    if "#" in i:
        pass
    else:
        vcf = re.split('[\t\n;]',i)
        x = str(vcf[7])
        x = str(vcf[-2]) 
        for c in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_transcript_ID.txt"):
            IDs = re.split('[\t\n]',c)
            if str(IDs[0]) in x: 
                transcript = str(vcf[4]) + "|" + str(IDs[0])
                scores = re.split('[\t\n,]',x)
                for score in scores:
                    if str(transcript) in score:
                        SpliceAI = re.split('[\t\n|]',score)
                        value = float(SpliceAI[4])
                        if float(value) >= 0.50:
                            value = "{0:.2f}".format(float(value))
                            outfile.writelines(i[:-1] + "\t" + str(value) +  "\t" + str(IDs[0]) + "\n")
                        if float(value) >= 0.80:
                            value = "{0:.2f}".format(float(value))
                            outfile2.writelines(i[:-1] + "\t" + str(value) +  "\t" + str(IDs[0]) + "\n")

outfile.close()
outfile2.close()

#------------------
#Donor site (T)
#------------------
print ("donoT")
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf","w")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over05.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai.vcf"):
    if "#" in i:
        pass
    else:
        vcf = re.split('[\t\n;]',i)
        x = str(vcf[7])
        x = str(vcf[-2]) 
        for c in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_transcript_ID.txt"):
            IDs = re.split('[\t\n]',c)
            if str(IDs[0]) in x: 
                transcript = str(vcf[4]) + "|" + str(IDs[0])
                scores = re.split('[\t\n,]',x)
                for score in scores:
                    if str(transcript) in score:
                        SpliceAI = re.split('[\t\n|]',score)
                        value = float(SpliceAI[4])
                        if float(value) >= 0.50:
                            value = "{0:.2f}".format(float(value))
                            outfile.writelines(i[:-1] + "\t" + str(value) +  "\t" + str(IDs[0]) + "\n")
                        if float(value) >= 0.80:
                            value = "{0:.2f}".format(float(value))
                            outfile2.writelines(i[:-1] + "\t" + str(value) +  "\t" + str(IDs[0]) + "\n")

outfile.close()
outfile2.close()


cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over05.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over05.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over05.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over05.vcf > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_omimgene_spliceai_exonslop50_over05.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)


cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_omimgene_spliceai_exonslop50_over08.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)
