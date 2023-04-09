#!/usr/bin/python

import sys
import re
import subprocess
import os
import gzip




outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/gnomad.genomes.r3.0.sites.all_oCDS_P3_allcat.getfasta_v2_allchromesome_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_PSI_domain_gene_refresh.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/gnomad.genomes.r3.0.sites.all_oCDS_P3_allcat.getfasta_v2_allchromesome_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_PSI_domain_gene.vcf"):
	vcf =  re.split('[\t\n]',i)
	isNothing = True
	for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/gnomad.genomes.r3.0.sites.allchr_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon_num_kai_erroruniq.vcf"):
		error = re.split('[\t\n]',c)
		if str(vcf[0]) == str(error[1]) and str(vcf[1]) == str(error[2]) and str(vcf[3]) == str(error[4]) and str(vcf[4]) == str(error[5]):
			isNothing = False
			if "#" in str(error[0]):
				pass
			if "#" not in str(error[0]) :
				outfile.writelines(str(vcf[0]) + "\t" + str(vcf[1]) + "\t" + str(vcf[2]) + "\t" + str(vcf[3]) + "\t" + str(vcf[4]) + "\t" + str(vcf[5]) + "\t" + str(vcf[6]) + "\t" + str(vcf[7]) + "\t" + str(vcf[8]) + "\t" + "HLA-C" + "\t" + str(error[10]) + "\t" + str(vcf[11]) + "\t" + str(vcf[12]) + "\t" + str(error[13]) + "\t" + str(error[14]) + "\t" + str(vcf[15]) + "\t" +  str(vcf[16]) + "\t" + str(vcf[17]) + "\t" + str(vcf[18]) + "\n")
				#print (i)
	if isNothing:
		outfile.writelines(i)

outfile.close()
