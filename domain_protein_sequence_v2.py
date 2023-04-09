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




#-------------------
#No PTC or Do not trigger NMD dirction
#-------------------



outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon_num_kai.vcf"):
	scm = re.split('[\t\n]',i)
	isNothing = True
	if "No PTC" in str(scm[14])  or "Do not trigger NMD" in str(scm[14]):
		#isNothing = False
		for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change.gtf"):
			gtf = re.split('[\t\n]',c)
			outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(scm[13]) + "\t" + str(scm[14]) + "\t" + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" + str(scm[18]) + "\t" + str(scm[19]) + "\t" + str(gtf[6]) + "\n")
			isNothing = False
			break
	if isNothing:
		outfile.writelines(i)


outfile.close()


#-------------------------------
#Get SCM gtf
#-------------------------------
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist.vcf"):
	scm = re.split('[\t\n]',i)
	if "No PTC" in str(scm[14])  and "pass" != str(scm[19]) : #"Do not trigger NMD" in str(scm[14]):
		outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change_codonpos.gtf","w")
		if "+" == str(scm[-2]):
			amari = 0
			sumlength = 0
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0])  + "_scm_change.gtf"):
				gtf = re.split('[\t\n]',d)
				transid = re.split('[\t\n;]',str(gtf[8]))
				if "CDS" == str(gtf[2]) :
					length  = int(gtf[4]) - int(gtf[3]) + 1 + amari
					amari = length % 3
					codon = length // 3
					sumlength +=  codon
					outfile.writelines(str(d[:-1]) + "\t" + str(transid[6]) + "\t" + str(length) + "\t" +  str(codon) + "\t" + str(sumlength) + "\n")

		if "-" == str(scm[-2]):
			amari = 0
			sumlength = 0
			for d in reversed(list(open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change.gtf"))):
				gtf = re.split('[\t\n]',d)
				transid = re.split('[\t\n;]',str(gtf[8]))
				if "CDS" == str(gtf[2]) :
					length  = int(gtf[4]) - int(gtf[3]) + 1 + amari
					amari = length % 3
					codon = length // 3
					sumlength += codon
					outfile.writelines(str(d[:-1]) + "\t" + str(transid[6]) + "\t" + str(length) + "\t" + str(codon) + "\t" + str(sumlength) + "\n")
		outfile.close()
	if "Do not trigger NMD" in str(scm[14])  and "pass" != str(scm[19]) :
		outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change_codonpos.gtf","w")
		if "+" == str(scm[-2]):
			amari = 0
			sumlength = 0
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0])  + "_scm_change.gtf"):
				gtf = re.split('[\t\n]',d)
				transid = re.split('[\t\n;]',str(gtf[8]))
				if "CDS" == str(gtf[2]) :
					length  = int(gtf[4]) - int(gtf[3]) + 1 + amari
					amari = length % 3
					codon = length // 3
					sumlength +=  codon
					outfile.writelines(str(d[:-1]) + "\t" + str(transid[6]) + "\t" + str(length) + "\t" +  str(codon) + "\t" + str(sumlength) + "\n")
		if "-" == str(scm[-2]):
			amari = 0
			sumlength = 0
			for d in reversed(list(open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change.gtf"))):
				gtf = re.split('[\t\n]',d)
				transid = re.split('[\t\n;]',str(gtf[8]))
				if "CDS" == str(gtf[2]) :
					length  = int(gtf[4]) - int(gtf[3]) + 1 + amari
					amari = length % 3
					codon = length // 3
					sumlength += codon
					outfile.writelines(str(d[:-1]) + "\t" + str(transid[6]) + "\t" + str(length) + "\t" + str(codon) + "\t" + str(sumlength) + "\n")
		outfile.close()
	







for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist.vcf"):
        scm = re.split('[\t\n]',c)
        if "No PTC" in str(scm[14])  or "Do not trigger NMD" in str(scm[14]) :
                if  "pass" != str(scm[19]):
                        #print (c)
                        hairetsu = str(scm[-3])
                        start = 0
                        outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change_codonpos_sequence.gtf","w")
                        for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change_codonpos.gtf"):
                                gtf = re.split('[\t\n]',d)
                                if "CDS" == str(gtf[2]):
                                        num = int(gtf[-3])
                                        end = start + num
                                        outfile.writelines(d[:-1] + "\t" + str(hairetsu[start:end]) + "\n")
                                        start = end
                        outfile.close()



#----------------------
#Annotation gtf
#----------------------




for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist.vcf"):
	scm = re.split('[\t\n]',c)
	if "No PTC" in str(scm[14])  or "Do not trigger NMD" in str(scm[14]) :
		if "pass" != str(scm[19]):
			i#print (c)
			outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre.gtf","w")
			outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_CDS.gtf","w")
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre.gtf"):
				gtf = re.split('[\t\n]',d)
				trans = re.split('[\t\n;]',str(gtf[9]))
				transID =re.split('[\t\n"]',str(trans[1]))
				if str(scm[10]) == str(transID[1]) and str(scm[0]) == str(gtf[0]):
					outfile.writelines(str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" + str(gtf[4]) + "\t" + str(gtf[5]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + str(gtf[8]) + "\t" + str(gtf[9]) + "\n")
					if "CDS" == str(gtf[3]):
						outfile2.writelines(str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" + str(gtf[4]) + "\t" + str(gtf[5]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + str(gtf[8]) + "\t" + str(gtf[9]) + "\n")
			outfile.close()
			outfile2.close()



#-------------
#gffread
#-------------


for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist.vcf"):
	scm = re.split('[\t\n]',c)
	if "No PTC" in str(scm[14])  or "Do not trigger NMD" in str(scm[14]) :
		if "pass" != str(scm[19]):
			print (c)
			cmd = "gffread -E /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre.gtf -g /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa -o /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_gffread.gff3"
			print (cmd)
			contents = subprocess.check_call(cmd,shell=True)

			cmd = "gffread /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_gffread.gff3 -g /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa -w /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_gffread.transcripts.fa -x /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_gffread.cds.fa -y /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_gffread.protein.fa"
			print (cmd)
			contents = subprocess.check_call(cmd,shell=True)



#-------------------------------
#Get annnotation gtf
#-------------------------------
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist.vcf"):
        scm = re.split('[\t\n]',i)
        if "No PTC" in str(scm[14])  and "pass" != str(scm[19]) : #"Do not trigger NMD" in str(scm[14]):
                outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_CDS_codonpos.gtf","w")
                if "+" == str(scm[-2]):
                        amari = 0
                        sumlength = 0
                        for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0])  + "_scm_pre_CDS.gtf"):
                                gtf = re.split('[\t\n]',d)
                                transid = re.split('[\t\n;]',str(gtf[8]))
                                if "CDS" == str(gtf[2]) :
                                        length  = int(gtf[4]) - int(gtf[3]) + 1 + amari
                                        amari = length % 3
                                        codon = length // 3
                                        sumlength +=  codon
                                        outfile.writelines(str(d[:-1]) + "\t" + str(transid[6]) + "\t" + str(length) + "\t" +  str(codon) + "\t" + str(sumlength) + "\n")

                if "-" == str(scm[-2]):
                        amari = 0
                        sumlength = 0
                        for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_CDS.gtf"):
                                gtf = re.split('[\t\n]',d)
                                transid = re.split('[\t\n;]',str(gtf[8]))
                                if "CDS" == str(gtf[2]) :
                                        length  = int(gtf[4]) - int(gtf[3]) + 1 + amari
                                        amari = length % 3
                                        codon = length // 3
                                        sumlength += codon
                                        outfile.writelines(str(d[:-1]) + "\t" + str(transid[6]) + "\t" + str(length) + "\t" + str(codon) + "\t" + str(sumlength) + "\n")
                outfile.close()
        if "Do not trigger NMD" in str(scm[14])  and "pass" != str(scm[19]) :
                outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_CDS_codonpos.gtf","w")
                if "+" == str(scm[-2]):
                        amari = 0
                        sumlength = 0
                        for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0])  + "_scm_pre_CDS.gtf"):
                                gtf = re.split('[\t\n]',d)
                                transid = re.split('[\t\n;]',str(gtf[8]))
                                if "CDS" == str(gtf[2]) :
                                        length  = int(gtf[4]) - int(gtf[3]) + 1 + amari
                                        amari = length % 3
                                        codon = length // 3
                                        sumlength +=  codon
                                        outfile.writelines(str(d[:-1]) + "\t" + str(transid[6]) + "\t" + str(length) + "\t" +  str(codon) + "\t" + str(sumlength) + "\n")
                if "-" == str(scm[-2]):
                        amari = 0
                        sumlength = 0
                        for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_CDS.gtf"):
                                gtf = re.split('[\t\n]',d)
                                transid = re.split('[\t\n;]',str(gtf[8]))
                                if "CDS" == str(gtf[2]) :
                                        length  = int(gtf[4]) - int(gtf[3]) + 1 + amari
                                        amari = length % 3
                                        codon = length // 3
                                        sumlength += codon
                                        outfile.writelines(str(d[:-1]) + "\t" + str(transid[6]) + "\t" + str(length) + "\t" + str(codon) + "\t" + str(sumlength) + "\n")
                outfile.close()
        




outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist_NO_PTC.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist.vcf"):
        scm = re.split('[\t\n]',c)
        if "No PTC" in str(scm[14])  or "Do not trigger NMD" in str(scm[14]) :
                if  "pass" != str(scm[19]):
                        hairetsu = ""
                        for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_gffread.protein.fa"):
                               if ">" in d:
                                       pass
                               else:
                                       amino = re.split('[\t\n]',d)
                                       hairetsu +=str(amino[0])
                        outfile.writelines(c[:-1] + "\t" + str(hairetsu) +"\n")




outfile.close()



for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist_NO_PTC.vcf"):
        scm = re.split('[\t\n]',c)
        if "No PTC" in str(scm[14])  or "Do not trigger NMD" in str(scm[14]) :
                if  "pass" != str(scm[19]):
                        hairetsu = str(scm[-2])
                        start = 0
                        outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_CDS_codonpos_sequence.gtf","w")
                        for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_CDS_codonpos.gtf"):
                                gtf = re.split('[\t\n]',d)
                                if "CDS" == str(gtf[2]):
                                        num = int(gtf[-3])
                                        end = start + num
                                        outfile.writelines(d[:-1] + "\t" + str(hairetsu[start:end]) + "\n")
                                        start = end
                        outfile.close()


outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist_NO_PTC_diffseq.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist_NO_PTC.vcf"):
	scm = re.split('[\t\n]',c)
	exonnum = re.split('[\t\n\s]',str(scm[15]))
	print (c)
	if "No PTC" in str(scm[14])  or "Do not trigger NMD" in str(scm[14]) :
		if  "pass" != str(scm[19]):
			pre = []
			aft = []
			outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_pre.diff_proteins.fasta","w")
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_pre_CDS_codonpos_sequence.gtf"):
				gtf = re.split('[\t\n]',d)
				gtfnum = re.split('[\t\n\s]',str(gtf[9]))
				if int(exonnum[2]) <= int(gtfnum[2]):
					pre.append(str(gtf[9]) + ":" + str(gtf[13]))
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change_codonpos_sequence.gtf"):
				gtf = re.split('[\t\n]',d)
				gtfnum = re.split('[\t\n\s]',str(gtf[9]))
				if int(exonnum[2]) <= int(gtfnum[2]):
					aft.append(str(gtf[9]) + ":" + str(gtf[13]))
			#print (aft)
			#print (len(aft))
			#if len(pre) == 1 and len(aft) == 1:
			diff_pre = ""
			diff_aft = ""
			pre_allsequence = ""
			aft_allsequence = ""
			for preseq, aftseq in zip(pre, aft):
				if preseq != aftseq:
					#print (aft)
					#print (len(aft))
					prehairetu = re.split('[\t\n:]',preseq)
					afthairetu = re.split('[\t\n:]',aftseq)
					diff_pre += str(prehairetu[0]) + ","
					diff_aft += str(afthairetu[0]) + ","
					pre_allsequence += str(prehairetu[1])
					aft_allsequence += str(afthairetu[1])
					#print (prehairetu)
					#print (afthairetu)
			outfile.writelines(">" + str(scm[10]) + "\n")
			outfile.writelines(pre_allsequence + "\n")
			outfile2.writelines(c[:-1] + "\t" +  str(diff_pre[:-1]) + "\t" + str(pre_allsequence) + "\n")
			print (aft_allsequence)
			#print (len(diff_pre),len(diff_aft))
			outfile.close()

outfile2.close()



#------------------------------
#domain or "Inter-domain region
#------------------------------


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist_NO_PTC_domain.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist_NO_PTC_diffseq.vcf"):
	scm = re.split('[\t\n]',c)
	exonnum = re.split('[\t\n\s]',str(scm[15]))
	if "No PTC" in str(scm[14])  or "Do not trigger NMD" in str(scm[14]) :
		if  "pass" != str(scm[19]):
			cmd = "hmmscan --tblout /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_pre.diff_proteins_hmmscan.tblout -E 0.0001 /mnt/data6/narumi/Pfam/Pfam-A.hmm /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_pre.diff_proteins.fasta"
			print (cmd)
			contents = subprocess.check_call(cmd,shell=True)
			hako = ""
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_pre.diff_proteins_hmmscan.tblout"):
				if "#" not in d:
					hmmscan = re.split('[\t\n\s]',d)
					hako += str(hmmscan[0]) + ","
			if len(hako) > 0:
				outfile.writelines(c[:-1] + "\t" + str(hako[:-1]) + "\n")
			if len(hako) == 0:
				outfile.writelines(c[:-1] + "\t" + "Inter-domain region" + "\n")

outfile.close()



#-----------------------
#Add domain infomation
#----------------------

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist_domain.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist.vcf"):
        scm = re.split('[\t\n]',c)
        isNothing = True
        for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist_NO_PTC_domain.vcf"):
                domain = re.split('[\t\n]',d)
                if str(scm[0]) == str(domain[0]) and str(scm[1]) == str(domain[1]) and str(scm[2]) == str(domain[2]) and str(scm[3]) == str(domain[3]) and str(scm[4]) == str(domain[4]) and str(scm[5]) == str(domain[5]):
                        isNothing = False
                        outfile.writelines(c[:-1] + "\t" + str(domain[-2]) + "\n")

        if isNothing:
                outfile.writelines(c[:-1] + "\t" + "-" + "\n")

outfile.close()


#----------
#sort col 
#-----------


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_PSI_domain.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_PSI.vcf"):
	scm = re.split('[\t\n]',c)
	isNothing = True
	for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_dirlist_domain.vcf"):
		domain = re.split('[\t\n]',d)
		if str(scm[0]) == str(domain[0]) and str(scm[1]) == str(domain[1]) and str(scm[2]) == str(domain[2]) and str(scm[3]) == str(domain[3]) and str(scm[4]) == str(domain[4]) and str(scm[5]) == str(domain[5]):
			if "No PTC" == str(domain[14]) or "Do not trigger NMD" == str(domain[14]):
				outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(scm[13]) + "\t" + str(domain[14])  + "\t" + str(domain[-2])  + "\t" + str(scm[14]) + "\t" + str(scm[15]) + "\t" +  str(scm[16]) + "\n")
			if "Not in CDS" == str(domain[14]) :
				outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(scm[13]) + "\t" + str(domain[14])  + "\t" + str(domain[-2])  + "\t" + "-" + "\t" + "-" + "\t" +  str(scm[16]) + "\n")

			if "No PTC" != str(domain[14]) and "Do not trigger NMD" != str(domain[14]) and  "Not in CDS" != str(domain[14]): 
				outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(scm[13]) + "\t" + str(domain[14])  + "\t" + str(domain[-2])  + "\t" + str(scm[14]) + "\t" + str(scm[15]) + "\t" +  str(scm[16]) + "\n")

outfile.close()



#--------------
#Add gene name
#--------------

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_PSI_domain_gene.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_GTEx_PSI_domain.vcf"):
	scm = re.split('[\t\n]',c)
	for d in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_transcript.gtf"):
		gtf = re.split('[\t\n]',d)
		detgtf = re.split('[\t\n;]',str(gtf[8]))
		transID = re.split('[\t\n"]',str(detgtf[1]))
		if str(transID[1]) == str(scm[10]) :
			genename =  re.split('[\t\n"]',str(detgtf[3]))
			outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(genename[1]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(scm[13]) + "\t" + str(scm[14])  + "\t" + str(scm[15])  + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" +  str(scm[18]) + "\n")
			break


outfile.close()


