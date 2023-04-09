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

import statistics
import math
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

'''
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_cp.vcf"):
	if "GTEX" in i:
		snv = re.split('[\t\n]',i)
		if "validated" == str(snv[14]):
			ids = re.split('[\t,]',str(snv[13]))
			for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon_num_kai.vcf"):
				hairetsu = re.split('[\t\n]',c)
				if str(snv[0]) == str(hairetsu[1]) and  str(snv[1]) == str(hairetsu[2]) and  str(snv[2]) == str(hairetsu[3]) and str(snv[3]) == str(hairetsu[4]) and str(snv[4]) == str(hairetsu[5]) and str(snv[5]) == str(hairetsu[6]):
					for IDs in ids:
						outfile.writelines(str(hairetsu[1]) + "\t" + str(hairetsu[16]) + "\t" + str(hairetsu[17]) + "\t" + str(hairetsu[10]) + "\t" + str(hairetsu[0]) + "\t" + str(hairetsu[20]) + "\t" + str(IDs) + "\t" + str(hairetsu[6]) + "\t" + str(hairetsu[13]) + "\t" + str(hairetsu[15]) + "\n")

outfile.close()


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_dir.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon.bed"):
	isNothing = True
	if "pass" in i:
		snv = re.split('[\t\n]',i)
		for c in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_" + chr + ".gtf"):
			if str(snv[3]) in c :
				isNothing = False
				gtf = re.split('[\t\n]',c)
				outfile.writelines(str(snv[0]) + "\t" + str(snv[1]) + "\t" + str(snv[2]) + "\t" + str(snv[3]) + "\t" + str(snv[4]) + "\t" + str(gtf[6]) + "\t" +  str(snv[6]) + "\t" +  str(snv[7]) + "\t" + str(snv[8]) + "\t" + str(snv[9]) + "\n")
				break
	if isNothing:
		outfile.writelines(i)


outfile.close()



for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_dir.bed"):
	snv = re.split('[\t\n]',c)
	mojiretsu =''
	mojiretsu = "$3~/" + str(snv[0]) + "/"
	cmd = "samtools view /mnt/data6/narumi/GTEx_analysis/blood/" + str(snv[6]) + "/hg38/hisat2/" + str(snv[6]) + ".hisat2.sorted.bam|awk '{if(" + mojiretsu + ")  print}' > /mnt/data6/narumi/GTEx_analysis/blood/" + str(snv[6]) + "/hg38/hisat2/" + str(snv[6]) + ".hisat2.scm_" + chr + ".sorted.bam"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)

	outfile = open("/mnt/data6/narumi/GTEx_analysis/blood/" + str(snv[6]) + "/hg38/hisat2/" + str(snv[6]) + ".hisat2.scm_" + chr + ".sorted2.bam","w")
	for i in open("/mnt/data6/narumi/GTEx_analysis/blood/" + str(snv[6]) + "/hg38/hisat2/" + str(snv[6]) + ".hisat2.scm_" + chr + ".sorted.bam"):
		bam =re.split('[\t\n]',i)
		if "N" in bam[5] and "S" not in bam[5] and "I" not in bam[5] and "D" not in bam[5] :
			outfile.writelines(i)

	outfile.close()



outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_junction.bam","w")

for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_dir.bed"):
	snv = re.split('[\t\n]',c)
	#print (c)
	for i in open("/mnt/data6/narumi/GTEx_analysis/blood/" + str(snv[6]) + "/hg38/hisat2/" + str(snv[6]) + ".hisat2.scm_" + chr + ".sorted2.bam"):
		situa =re.split('[\t\n]',i)
		dd = format(int(situa[1]),'b')
		cigar = re.split('(M|N)',situa[5])
		if "M" in cigar[1] and "N" in cigar[3] and str(snv[0]) == str(situa[2]):
			if dd[-5] in "1" or dd[-5] in "0":
				if "+" in str(snv[5]) and int(snv[1]) == int(situa[3]) + int (cigar[0]) + int (cigar[2]) and "3′ss" == str(snv[7]):
					outfile.writelines(str(situa[0]) +"\t" + str(situa[1]) +"\t" + str(situa[2]) +"\t" +str(situa[3]) +"\t" + str(situa[4]) +"\t" + str(situa[5]) +"\t" + str(situa[6]) +"\t" + str(situa[7]) +"\t" + str(situa[8]) + "\t" + "3′ss" + "\t" + "novel exon" + "\t" + str(snv[8]) + "\t" + str(snv[4]) + "\t" + str(snv[6]) +"\n")
				if "+" in snv[5] and int(snv[2]) == int(situa[3]) + int (cigar[0])  -1 and "5′ss" == str(snv[7]):
					outfile.writelines(str(situa[0]) +"\t" + str(situa[1]) +"\t" + str(situa[2]) +"\t" +str(situa[3]) +"\t" + str(situa[4]) +"\t" + str(situa[5]) +"\t" + str(situa[6]) +"\t" + str(situa[7]) +"\t" + str(situa[8]) + "\t" + "5′ss" + "\t" + "novel exon" + "\t"+ str(snv[8]) + "\t" + str(snv[4]) + "\t" + str(snv[6]) + "\n")
				if "-" in str(snv[5]) and int(snv[2]) == int(situa[3]) + int (cigar[0]) - 1 and "3′ss" == str(snv[7]):
					outfile.writelines(str(situa[0]) +"\t" + str(situa[1]) +"\t" + str(situa[2]) +"\t" +str(situa[3]) +"\t" + str(situa[4]) +"\t" + str(situa[5]) +"\t" + str(situa[6]) +"\t" + str(situa[7]) +"\t" + str(situa[8]) + "\t" + "3′ss" + "\t" + "novel exon" + "\t" + str(snv[8]) + "\t" + str(snv[4]) + "\t" + str(snv[6]) + "\n")
				if "-" in str(snv[5]) and int(snv[1]) == int(situa[3]) + int (cigar[0])  +  int (cigar[2]) and "5′ss" == str(snv[7]):
					outfile.writelines(str(situa[0]) +"\t" + str(situa[1]) +"\t" + str(situa[2]) +"\t" +str(situa[3]) +"\t" + str(situa[4]) +"\t" + str(situa[5]) +"\t" + str(situa[6]) +"\t" + str(situa[7]) +"\t" + str(situa[8]) + "\t" + "5′ss" + "\t" + "novel exon" + "\t" + str(snv[8]) + "\t" + str(snv[4]) + "\t" + str(snv[6]) + "\n")

outfile.close()

					

#-----------------------
#annotation splice site
#-----------------------

outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.annotationexon.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_dir.bed"):
	snv = re.split('[\t\n]',i)
	for c in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_" + chr + ".gtf"):
		if str(snv[3]) in c and str(snv[9]) in c:
			gtf = re.split('[\t\n]',c)
			outfile2.writelines(str(gtf[0]) + "\t" + str(gtf[3]) + "\t" + str(gtf[4]) + "\t" + str(snv[3]) + "\t" + str(snv[4]) +  "\t" + str(gtf[6]) + "\t" +  str(snv[6]) + "\t" +  str(snv[7]) + "\t" + str(snv[8]) + "\t" + str(snv[9]) + "\n")
			break


outfile2.close()


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.annotationexon_junction.bam","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.annotationexon.bed"):
	snv = re.split('[\t\n]',c)
	#print (c)
	for i in open("/mnt/data6/narumi/GTEx_analysis/blood/" + str(snv[6]) + "/hg38/hisat2/" + str(snv[6]) + ".hisat2.scm_" + chr + ".sorted2.bam"):
		situa =re.split('[\t\n]',i)
		dd = format(int(situa[1]),'b')
		cigar = re.split('(M|N)',situa[5])
		if "M" in cigar[1] and "N" in cigar[3] and str(snv[0]) == str(situa[2]):
			if dd[-5] in "1" or dd[-5] in "0":
				if "+" in str(snv[5]) and int(snv[1]) == int(situa[3]) + int (cigar[0]) + int (cigar[2]) and "3′ss" == str(snv[7]):
					outfile.writelines(str(situa[0]) +"\t" + str(situa[1]) +"\t" + str(situa[2]) +"\t" +str(situa[3]) +"\t" + str(situa[4]) +"\t" + str(situa[5]) +"\t" + str(situa[6]) +"\t" + str(situa[7]) +"\t" + str(situa[8]) + "\t" + "3′ss" + "\t" + "annotation exon" + "\t" + str(snv[8]) + "\t" + str(snv[4]) + "\t" + str(snv[6])  +"\n")
				if "+" in snv[5] and int(snv[2]) == int(situa[3]) + int (cigar[0])  -1 and "5′ss" == str(snv[7]):
					outfile.writelines(str(situa[0]) +"\t" + str(situa[1]) +"\t" + str(situa[2]) +"\t" +str(situa[3]) +"\t" + str(situa[4]) +"\t" + str(situa[5]) +"\t" + str(situa[6]) +"\t" + str(situa[7]) +"\t" + str(situa[8]) + "\t" + "5′ss" + "\t" + "annotation exon" + "\t" + str(snv[8]) + "\t" + str(snv[4]) + "\t" + str(snv[6]) + "\n")
				if "-" in str(snv[5]) and int(snv[2]) == int(situa[3]) + int (cigar[0]) - 1 and "3′ss" == str(snv[7]):
					outfile.writelines(str(situa[0]) +"\t" + str(situa[1]) +"\t" + str(situa[2]) +"\t" +str(situa[3]) +"\t" + str(situa[4]) +"\t" + str(situa[5]) +"\t" + str(situa[6]) +"\t" + str(situa[7]) +"\t" + str(situa[8]) + "\t" + "3′ss" + "\t" + "annotation exon" + "\t" + str(snv[8]) + "\t" + str(snv[4]) + "\t" + str(snv[6]) + "\n")
				if "-" in str(snv[5]) and int(snv[1]) == int(situa[3]) + int (cigar[0])  +  int (cigar[2]) and "5′ss" == str(snv[7]):
					outfile.writelines(str(situa[0]) +"\t" + str(situa[1]) +"\t" + str(situa[2]) +"\t" +str(situa[3]) +"\t" + str(situa[4]) +"\t" + str(situa[5]) +"\t" + str(situa[6]) +"\t" + str(situa[7]) +"\t" + str(situa[8]) + "\t" + "5′ss" + "\t" + "annotation exon" + "\t" + str(snv[8]) + "\t" + str(snv[4]) + "\t" + str(snv[6]) + "\n")
outfile.close()





#-------------------
#Caculated PSI value
#--------------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_PSIvalue.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_dir.bed"):
	snv = re.split('[\t\n]',c)
	novelcount =[]
	#print (c)
	for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_junction.bam"):
		bam = re.split('[\t\n]',i)
		#print (snv)
		#print (bam)
		if str(snv[4]) == str(bam[12]) and str(snv[6]) == str(bam[13]) :
			novelcount.append(i)
	#print (len(novelcount))
	annocount =[]
	for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.annotationexon_junction.bam"):
		bam = re.split('[\t\n]',i)
		if str(snv[4]) == str(bam[12]) and str(snv[6]) == str(bam[13]) :
			annocount.append(i)
	#print (len(annocount))
	if "In intronic" == str(snv[8]) :
		PSI = float(len(novelcount)) / float(len(novelcount) + len(annocount))
		#if PSI > 0:
		outfile.writelines(str(c[:-1]) + "\t" + str(len(novelcount)) +  "\t" + str(len(annocount)) + "\t" + str(PSI) + "\n")
	if "In intronic" != str(snv[8]) and "3′UTR" != str(snv[8]) and "5′UTR" != str(snv[8]) and len(annocount) > 0:
		PSI = float(len(annocount)) / float(len(novelcount) + len(annocount))
		#if PSI > 0:
		outfile.writelines(str(c[:-1]) + "\t" + str(len(annocount)) +  "\t" + str(len(novelcount)) + "\t" + str(PSI) + "\n")
	if "In intronic" != str(snv[8]) and "3′UTR" != str(snv[8]) and "5′UTR" != str(snv[8]) and len(annocount) == 0:
		outfile.writelines(str(c[:-1]) + "\t" + str(len(annocount)) +  "\t" + str(len(novelcount)) + "\t" + "0.00" + "\n")
	if "3′UTR" == str(snv[8]) or "5′UTR" == str(snv[8]):
		outfile.writelines(str(c[:-1]) + "\t" + str(len(annocount)) +  "\t" + str(len(novelcount)) + "\t" + "Not caculated" + "\n")


outfile.close()
'''
#------------------------------
#PSI mean & standard deviation
#------------------------------

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_PSI0.vcf","w")
num = 0
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx.vcf"):
	num += 1
	snv = re.split('[\t\n]',c)
	IDs = ""
	psival = []
	if "3′UTR" != str(snv[12]) and "5′UTR" != str(snv[12]):
		for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_PSIvalue.bed"):
			psi = re.split('[\t\n]',d)
			if str(num) == str(psi[4]) :
				IDs += str(psi[6])+ ","
				psival.append(float(psi[12]))
		if len(psival) > 1:
			mean = statistics.mean(psival)
			pstdev = statistics.pstdev(psival)
			outfile.writelines(str(num) + "\t" + c[:-1] + "\t" + str(IDs[:-1]) + "\t" + str(round(mean,2)) + "±" + str(round(pstdev,2)) + "\n")
		if len(psival) == 1:
			mean = statistics.mean(psival)
			outfile.writelines(str(num) + "\t" + c[:-1] + "\t" + str(IDs[:-1]) + "\t" + str(round(mean,2)) + "\n")
		if len(psival) == 0 :
			if "GTEX-" in str(snv[13]):
				outfile.writelines(str(num) + "\t" + c[:-1] + "\t" + "-" + "\t" + "-" + "\n")
			else:
				outfile.writelines(str(num) + "\t" + c[:-1] + "-" + "\t" + "-" + "\t" + "-" + "\n")
		#outfile.writelines(str(num) + "\t" + c[:-1] + "\t" + str(IDs[:-1]) + "\t" + "Not caculated"+ "\n")
	
	if "3′UTR" == str(snv[12]) or "5′UTR" == str(snv[12]):		
	#IDs = ""
		for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/" + chr + "_GTEx_Analysis.novelexon_PSIvalue.bed"):
			psi = re.split('[\t\n]',d)
			if str(num) == str(psi[4]): 
		#	if  "3′UTR" == str(psi[8]) or "5′UTR" == str(psi[8]) :
				IDs += str(psi[6])+ ","
		outfile.writelines(str(num) + "\t" + c[:-1] + "\t" + str(IDs[:-1]) + "\t" + "Not caculated"+ "\n")
	
			
outfile.close()

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_PSI.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_PSI0.vcf"):
	snv = re.split('[\t\n]',c)
	isNothing = True
	if "GTEX-" in str(snv[14]) and "-" == str(snv[15]):
		num = 0
		for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_cp.vcf"):
			num += 1
			cp = re.split('[\t\n]',d)
			if str(snv[0]) == str(num) and str(snv[1]) == str(cp[0]) and str(snv[2]) == str(cp[1])  and str(snv[3]) == str(cp[2]) and str(snv[4]) == str(cp[3]) and str(snv[5]) == str(cp[4]) :
				isNothing = False
				outfile.writelines(str(snv[0]) + "\t" + str(snv[1]) + "\t" + str(snv[2]) + "\t" + str(snv[3]) + "\t" + str(snv[4]) + "\t" + str(snv[5]) + "\t" + str(snv[6]) + "\t" + str(snv[7]) + "\t" + str(snv[8]) + "\t" + str(snv[9]) + "\t" + str(snv[10]) + "\t" + str(snv[11]) + "\t" + str(snv[12]) + "\t" + str(snv[13]) + "\t" + str(snv[14]) + "\t" + str(cp[14]) + "\t" + str(cp[15]) + "\n")

	if isNothing:
		outfile.writelines(c)

outfile.close()

