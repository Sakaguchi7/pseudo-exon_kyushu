#!/usr/bin/python

import sys
import re
import subprocess
import os
import gzip

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


#----------------
#Acceptor site (A)
#-----------------

print ("accA")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_snv.vcf"):
    vcf =  re.split('[\t\n]',i)
    details = re.split('[\t\n;]',str(vcf[7])) 
    AFs = re.split('[\t\n=]', str(details[2]))
    if "AF" == str(AFs[0]):
        #print (AFs)
        score = float(AFs[1])
        if score <= 0.01 :
            outfile.writelines(i)

outfile.close()

#----------------
#Acceptor site (G)
#-----------------
print ("accG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_snv.vcf"):
    vcf =  re.split('[\t\n]',i)
    details = re.split('[\t\n;]',str(vcf[7])) 
    AFs = re.split('[\t\n=]', str(details[2]))
    if "AF" == str(AFs[0]):
        #print (AFs)
        score = float(AFs[1])
        if score <= 0.01 :
            outfile.writelines(i)

outfile.close()


#----------------
#Donor site (G)
#-----------------
print ("donoG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_snv.vcf"):
    vcf =  re.split('[\t\n]',i)
    details = re.split('[\t\n;]',str(vcf[7])) 
    AFs = re.split('[\t\n=]', str(details[2]))
    if "AF" == str(AFs[0]):
        #print (AFs)
        score = float(AFs[1])
        if score <= 0.01 :
            outfile.writelines(i)

outfile.close()


#----------------
#Donor site (T)
#-----------------
print ("donoT")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_snv.vcf"):
    vcf =  re.split('[\t\n]',i)
    details = re.split('[\t\n;]',str(vcf[7])) 
    AFs = re.split('[\t\n=]', str(details[2]))
    if "AF" == str(AFs[0]):
        #print (AFs)
        score = float(AFs[1])
        if score <= 0.01 :
            outfile.writelines(i)

outfile.close()




