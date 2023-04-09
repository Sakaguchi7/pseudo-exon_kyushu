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



num = 0
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.addexonlength.vcf"):
	scm = re.split('[\t\n]',i)
	#print (i)
	num += 1
	hairetsu = ""
	isNothing = True
	if "In CDS" == str(scm[13]):
		isNothing = False
		for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(num) + "_scm_change_gffread.protein.fa"):
			if ">" in c:
				pass
			else:
				amino = re.split('[\t\n]',c)
				hairetsu +=str(amino[0])
		ptc_pos = hairetsu.find('.')
		if int(ptc_pos) == -1:
			outfile.writelines(str(num) + "\t" + str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + "No PTC" + "\t" + str(scm[14]) + "\t" + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" + str(hairetsu) +"\n")
		else:
			outfile.writelines(str(num) + "\t" + str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(ptc_pos) + "\t" + str(scm[14]) + "\t" + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" + str(hairetsu) +"\n")
	if isNothing:
		outfile.writelines(str(num) + "\t" + str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(scm[13]) + "\t" + str(scm[14]) + "\t" + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\n")

outfile.close()




outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir.vcf","w")
#num = 0
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence.vcf"):
	scm = re.split('[\t\n]',i)
	#num += 1
	isNothing = True
	#print (scm)
	if "Not in CDS" != str(scm[14]) and "Do not trigger NMD" != str(scm[14]) and "No PTC" != str(scm[14]) :
		#isNothing = False
		for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change.gtf"):
			gtf = re.split('[\t\n]',d)
			transid = re.split('[\t\n;]',str(gtf[8]))
			ids =  re.split('[\t\n"]',str(transid[1]))
			if str(ids[1]) == str(scm[10]) :
				isNothing = False
				outfile.writelines(i[:-1] + "\t" + str(gtf[6]) +  "\n")
				break
	if isNothing:
		outfile.writelines(i[:-1] + "\t" + "pass" +  "\n")

outfile.close()





outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_frameshift_CDS.gtf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_exonlength_pre_countcodon_max.gtf"):
	codon = re.split('[\t\n]',i)
	count = -1
	#print (i)
	for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre.gtf"):
		gtf = re.split('[\t\n]',d)
		if "CDS" == str(gtf[3]):
			if int(codon[1]) == int(gtf[0]):
				count += 1
				num = codon[0][count]
				if str(num) == str(gtf[8]):
					pass
				else:
					outfile.writelines(d)
					break

outfile.close()





	


#---------------
#
#----------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir.vcf"):
	scm = re.split('[\t\n]',i)
        #num += 1
        #isNothing = True
	#print (i)
	isNothing = True
	if "Not in CDS" != str(scm[14]) and "Do not trigger NMD" != str(scm[14]) and "No PTC" != str(scm[14]) :
		for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_frameshift_CDS.gtf"):
			gtf = re.split('[\t\n]',d)
			if str(scm[0]) == str(gtf[0]) :
				isNothing = False
				outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" +  str(scm[13]) + "\t" + "PTC by frameshift" + "\t"  + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" +  str(scm[18]) + "\t" + str(scm[19]) + "\t" +str(scm[20]) + "\n")
				break		

	if isNothing:
		outfile.writelines(i)

outfile.close()

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift.vcf"):
	scm = re.split('[\t\n]',i)
	isNothing = True
	if "Not in CDS" != str(scm[14]) and "Do not trigger NMD" != str(scm[14]) and "No PTC" != str(scm[14]) and "PTC by frameshift" != str(scm[14]):
			isNothing = False
			outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" +  str(scm[13]) + "\t" + "PTC in novel exon" + "\t"  + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" +  str(scm[18]) + "\t" + str(scm[19]) + "\t" +str(scm[20]) + "\n")

	if isNothing:
		outfile.writelines(i)
outfile.close()
	







#-------------------------------------
#PTC in novel exon or PTC by frameshift
#---------------------------------------


for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon.vcf"):
	scm = re.split('[\t\n]',i)
	if "PTC by frameshift" == str(scm[14]) or "PTC in novel exon" == str(scm[14]) :
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




outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon_num.vcf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon.vcf"):
	scm = re.split('[\t\n]',i)
	isNothing = True
	#print (scm)
	if "PTC by frameshift" == str(scm[14]) or "PTC in novel exon" == str(scm[14]) :
		exon = re.split('[\s\n]',str(scm[15]))
		ptc_pos = str(scm[-3]).find('.')
		codon_pos = int(ptc_pos) + 1
		#print (exon)
		#if "+" == str(scm[-2]) and "3′ss" == str(scm[6]) :
		if "-" == str(scm[-2]) or  "+" == str(scm[-2]) :
			for  d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0])  + "_scm_change_codonpos.gtf"):
				gtf = re.split('[\t\n]',d)
				det = re.split('[\t\n;]',str(gtf[8]))
				gtf_exon = re.split('[\s\n]',str(det[6]))
				if int(codon_pos) <= int(gtf[-2]):
					#if int(exon[2]) < int(gtf_exon[2]):
					#	break
					#print (str(exon),str(gtf_exon))
					if int(exon[2])	== int(gtf_exon[2]):
						isNothing = False
						outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(scm[13]) + "\t" + str(scm[14])  + "\t" + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" +  str(scm[18]) + "\t" + str(scm[19]) + "\t" +str(scm[20]) + "\t" + str(codon_pos) + "\t" + str(det[6]) + "\n")
						break
					if int(exon[2]) < int(gtf_exon[2]):
						isNothing = False
						outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" +  str(scm[13]) + "\t" + str(scm[14]) + "\t"  + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" +  str(scm[18]) + "\t" + str(scm[19]) + "\t" +str(scm[20]) + "\t" + str(codon_pos) + "\t" + str(det[6])  + "\n")
						break
					#if int(exon[2]) > int(gtf_exon[2]):
					#	print (str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" +  str(scm[13]) + "\t" + str(scm[14]) + "\t"  + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" +  str(scm[18]) + "\t" + str(scm[19]) + "\t" +str(scm[20]) + "\t" + str(codon_pos) + "\t" + str(det[6])  + "\n")
					#	isNothing = False
					#	outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" +  str(scm[13]) + "\t" + str(scm[14]) + "\t"  + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" +  str(scm[18]) + "\t" + str(scm[19]) + "\t" +str(scm[20]) + "\t" + str(codon_pos) + "\t" + str(det[6])  + "\n")
					#	break
		
	if isNothing:
		outfile.writelines(i)

outfile.close()




for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon_num.vcf"):
	scm = re.split('[\t\n]',c)
	if "PTC by frameshift" == str(scm[14]) or "PTC in novel exon" == str(scm[14]) :
		hairetsu = str(scm[-5])
		#print (hairetsu)
		start = 0
		outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change_codonpos_sequence.gtf","w")
		#for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_change.gtf"):
		for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change_codonpos.gtf"):
			gtf = re.split('[\t\n]',d)
			if "CDS" == str(gtf[2]):
				num = int(gtf[-3])
				end = start + num
				outfile.writelines(d[:-1] + "\t" + str(hairetsu[start:end]) + "\n")
				start = end
		outfile.close()
	
				


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon_num_kai.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon_num.vcf"):
	scm = re.split('[\t\n]',c)
	#print (scm)
	isNothing = True
	if "PTC by frameshift" == str(scm[14]) or "PTC in novel exon" == str(scm[14]) :
		isNothing = False
		for d in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_lastcds_number.txt"):
			numcds = re.split('[\t\n]',d)
			if str(numcds[0]) == str(scm[10]):
				last_exon =  re.split('[\t\s]', numcds[1])
				ptc_exon = re.split('[\t\s]', str(scm[-2]))
				#print (ptc_exon,last_exon)
				if int(last_exon[2]) == int(ptc_exon[2]):
					outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(scm[13]) + "\t" + "Do not trigger NMD" +  "\t" + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" +  str(scm[18]) + "\t" + str(scm[19]) + "\t" +str(scm[20]) + "\t" +  str(scm[21]) + "\t" +str(scm[22]) + "\n") 

				if int(last_exon[2]) - 1 == int(ptc_exon[2]):
					for f in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(scm[0]) + "_scm_change_codonpos_sequence.gtf"):
						gtf = re.split('[\t\n]',f)
						if str(gtf[9]) == str(scm[-2]):
							ptcpos = str(gtf[12]).find(".")
							if len(str(gtf[12])) - ptcpos < 16 :
								outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + str(scm[12]) + "\t" + str(scm[13]) + "\t" + "Do not trigger NMD" +  "\t" + str(scm[15]) + "\t" + str(scm[16]) + "\t" + str(scm[17]) + "\t" +  str(scm[18]) + "\t" + str(scm[19]) + "\t" +str(scm[20]) + "\t" +  str(scm[21]) + "\t" +str(scm[22]) + "\n") 
							else:
								outfile.writelines(c)
				if int(last_exon[2]) != int(ptc_exon[2]) and int(last_exon[2]) - 1 != int(ptc_exon[2]):
					outfile.writelines(c)
	if isNothing:
		outfile.writelines(c)


	
outfile.close()






#cmd = "less -S /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon_num_kai.vcf|cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,21,22,23 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.sequence_dir_frameshift_novelexon_num_kaicut.vcf"
#print (cmd)
#contents = subprocess.check_call(cmd,shell=True)

