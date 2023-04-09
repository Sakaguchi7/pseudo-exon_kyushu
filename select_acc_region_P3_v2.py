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

#-----------------------
#Acceptor site (A)
#------------------------

print ("accA")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.vcf"):
	if "#" in i:
		pass
	else:
		snp = re.split('[\t\n]',i)
		if "+" == str(snp[-2]):
			start = int(snp[1]) - 19
			end = int(snp[1]) + 4
			outfile.writelines(str(snp[0]) + "\t" + str(start) + "\t" + str(end) + "\t" + str(snp[3]) + "\t" + str(snp[4]) + "\t" + str(snp[-2]) + "\n")
		if "-" == str(snp[-2]):
			start = int(snp[1])	- 5
			end = int(snp[1]) + 18
			outfile.writelines(str(snp[0]) + "\t" + str(start) + "\t" + str(end) + "\t" + str(snp[3]) + "\t" + str(snp[4]) + "\t" + str(snp[-2]) + "\n")
outfile.close()

cmd = "bedtools getfasta -fi /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa -bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.bed -s -name -bedOut > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta.bed"):
	snp = re.split('[\t\n]',i)
	if "N" in str(snp[6]):
		pass
	else:
		outfile.writelines(i)
outfile.close()






#-----------------------
#Acceptor site (G)
#------------------------

print ("accG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoaccG.vcf"):
	if "#" in i:
		pass
	else:
		snp = re.split('[\t\n]',i)
		if "+" == str(snp[-2]):
			start = int(snp[1]) - 20
			end = int(snp[1]) + 3
			outfile.writelines(str(snp[0]) + "\t" + str(start) + "\t" + str(end) + "\t" + str(snp[3]) + "\t" + str(snp[4]) + "\t" + str(snp[-2]) +  "\n")
		if "-" == str(snp[-2]):
			start = int(snp[1]) - 4
			end = int(snp[1]) + 19
			outfile.writelines(str(snp[0]) + "\t" + str(start) + "\t" + str(end) + "\t" + str(snp[3]) + "\t" + str(snp[4]) + "\t" + str(snp[-2]) + "\n")
outfile.close()

cmd = "bedtools getfasta -fi /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa -bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.bed -s -name -bedOut > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta.bed"):
	snp = re.split('[\t\n]',i)
	if "N" in str(snp[6]):
		pass
	else:
		outfile.writelines(i)
outfile.close()


