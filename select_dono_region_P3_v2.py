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
#Donor site (T)
#------------------------
print (chr)
print ("donoT")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.vcf"):
	if "#" in i:
		pass
	else:
		snp = re.split('[\t\n]',i)
		if "+" == str(snp[-2]):
			start = int(snp[1]) - 5
			end = int(snp[1]) + 4
			outfile.writelines(str(snp[0]) + "\t" + str(start) + "\t" + str(end) + "\t" + str(snp[3]) + "\t" + str(snp[4]) + "\t" + str(snp[-2]) + "\n")
		if "-" == str(snp[-2]):
			start = int(snp[1])	- 5
			end = int(snp[1]) + 4
			outfile.writelines(str(snp[0]) + "\t" + str(start) + "\t" + str(end) + "\t" + str(snp[3]) + "\t" + str(snp[4]) + "\t" + str(snp[-2])  + "\n")
outfile.close()

cmd = "bedtools getfasta -fi /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa -bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.bed -s -name -bedOut > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta.bed"):
	snp = re.split('[\t\n]',i)
	if "N" in str(snp[6]):
		pass
	else:
		outfile.writelines(i)
outfile.close()




#-----------------------
#Donor site (G)
#------------------------

print ("donoG")	
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoaccG.vcf"):
	if "#" in i:
		pass
	else:
		snp = re.split('[\t\n]',i)
		if "+" == str(snp[-2]):
			start = int(snp[1]) - 4
			end = int(snp[1]) + 5
			outfile.writelines(str(snp[0]) + "\t" + str(start) + "\t" + str(end) + "\t" + str(snp[3]) + "\t" + str(snp[4]) + "\t" + str(snp[-2]) + "\n")
		if "-" == str(snp[-2]):
			start = int(snp[1])	- 6
			end = int(snp[1]) + 3
			outfile.writelines(str(snp[0]) + "\t" + str(start) + "\t" + str(end) + "\t" + str(snp[3]) + "\t" + str(snp[4]) + "\t" + str(snp[-2]) + "\n")
outfile.close()


cmd = "bedtools getfasta -fi /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa -bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.bed -s -name -bedOut > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta.bed"):
	snp = re.split('[\t\n]',i)
	if "N" in str(snp[6]):
		pass
	else:
		outfile.writelines(i)
outfile.close()


