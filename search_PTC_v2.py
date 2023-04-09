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


'''

#-----------------------------
#non-coding RNA or coding RNA
#------------------------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_lincRNA.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR.vcf"):
	scm =  re.split('[\t\n]',c)
	isNothing = True
	for d in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_transcript_ID_ncRNAlist.txt"):
		ncRNA = re.split('[\t\n]',d)
		if str(scm[9]) == str(ncRNA[0]):
			isNothing = False
			if "In intronic" == str(scm[12]):
				#isNothing = False
				outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + "lincRNA in intronic" + "\n")
			else:
				outfile.writelines(str(scm[0]) + "\t" + str(scm[1]) + "\t" + str(scm[2]) + "\t" + str(scm[3]) + "\t" + str(scm[4]) + "\t" + str(scm[5]) + "\t" + str(scm[6]) + "\t" + str(scm[7]) + "\t" + str(scm[8]) + "\t" + str(scm[9]) + "\t" + str(scm[10]) + "\t" + str(scm[11]) + "\t" + "lincRNA" + "\n")
	if isNothing:
		outfile.writelines(c)

outfile.close()


#-----------------------------
#Not in CDS or In CDS & lastEJC
#------------------------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.inCDS_EJC.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR_lincRNA.vcf"):
	scm = re.split('[\t\n]',c)
	isNothing = True
	#print (scm)
	for i in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_CDS_and_EJC_editv2.gff3"):
		cds = re.split('[\t\n]',i)
		#print (cds)
		if str(scm[9]) in str(cds[6]):
			isNothing = False
			#print (scm)
			if int(cds[3]) <= int(scm[1]) and int(cds[4]) >= int(scm[1]):
				if int(cds[9]) <= int(scm[1]) and int(cds[10]) >= int(scm[1]):
					outfile.writelines(c[:-1] + "\t" + "In CDS" + "\n")
					break
				if int(cds[9]) > int(scm[1]) or int(cds[10]) < int(scm[1]):
					outfile.writelines(c[:-1] + "\t" + "Do not trigger NMD" + "\n")
					#print (c[:-1] + "\t" + "Do not trigger NMD" + "\n")
					break
			if int(cds[3]) > int(scm[1]) or int(cds[4]) < int(scm[1]):
				outfile.writelines(c[:-1] + "\t" + "Not in CDS" + "\n")
				break
	if isNothing:
		outfile.writelines(c[:-1] + "\t" + "Not in CDS" + "\n")
outfile.close()


#-----------------
#Choose gtf data
#-----------------
cmd = "awk '{if($1==\"" + chr + "\"){print $0}}' /mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29.gtf > /mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_" + chr + ".gtf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

num = 0
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre.gtf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.inCDS_EJC.vcf"):
	scm = re.split('[\t\n]',c)
	num += 1
	if "In CDS" == str(scm[13]):
		for d in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_" + chr + ".gtf"):
			gtf = re.split('[\t\n]',d)
			trans = re.split('[\t\n;]',str(gtf[8]))
			transID =re.split('[\t\n"]',str(trans[1]))
			#print (transID)
			if str(scm[9]) == str(transID[1]):
				outfile.writelines(str(num) + "\t" + d)

outfile.close()


#------------
#Exonic SCM
#------------
cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_acceptorAG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)


cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donorGT.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

#-------------
#Intronic SCM
#-------------

cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_acceptorAG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donorGT.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)




#---------------
#Add exon number
#----------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.addexonlength.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.inCDS_EJC.vcf"):
	scm = re.split('[\t\n]',c)
	isNothing = True
	if "In intronic" != str(scm[12]):
		if "3′ss" == str(scm[5]):
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_acceptorAG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed"):
				newexon = re.split('[\t\n]',d)
				if str(scm[0]) == str(newexon[0]) and str(scm[1]) == str(newexon[4]) and str(scm[9]) == str(newexon[7]):
					isNothing = False
					outfile.writelines(c[:-1] + "\t" + str(newexon[3]) + "\t" + str(newexon[1]) + "\t" + str(newexon[2]) + "\t" + str(newexon[6]) + "\n")
					break
		if "5′ss" == str(scm[5]):
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donorGT.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed"):
				newexon = re.split('[\t\n]',d)
				if str(scm[0]) == str(newexon[0]) and str(scm[1]) == str(newexon[4]) and str(scm[9]) == str(newexon[7]):
					isNothing = False
					outfile.writelines(c[:-1] + "\t" + str(newexon[3]) + "\t" + str(newexon[1]) + "\t" + str(newexon[2]) + "\t" + str(newexon[6]) + "\n")
					break
	if "In intronic" == str(scm[12]):
		if "3′ss" == str(scm[5]):
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_acceptorAG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed"):
				newexon = re.split('[\t\n]',d)
				if str(scm[0]) == str(newexon[0]) and str(scm[1]) == str(newexon[5]) and str(scm[9]) == str(newexon[8]):
					isNothing = False
					outfile.writelines(c[:-1] + "\t" + str(newexon[4]) + "\t" + str(newexon[1]) + "\t" + str(newexon[2]) + "\t" + str(newexon[7]) + "\n")
					break
		if "5′ss" == str(scm[5]):
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donorGT.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed"):
				newexon = re.split('[\t\n]',d)
				if str(scm[0]) == str(newexon[0]) and str(scm[1]) == str(newexon[5]) and str(scm[9]) == str(newexon[8]):
					isNothing = False
					outfile.writelines(c[:-1] + "\t" + str(newexon[4]) + "\t" + str(newexon[1]) + "\t" + str(newexon[2]) + "\t" + str(newexon[7]) + "\n")
					break

	if "lincRNA in intronic" == str(scm[12]):
		if "3′ss" == str(scm[5]):
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_acceptorAG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed"):
				newexon = re.split('[\t\n]',d)
				if str(scm[0]) == str(newexon[0]) and str(scm[1]) == str(newexon[5]) and str(scm[9]) == str(newexon[8]):
					isNothing = False
					outfile.writelines(c[:-1] + "\t" + str(newexon[4]) + "\t" + str(newexon[1]) + "\t" + str(newexon[2]) + "\t" + str(newexon[7]) + "\n")
					break

		if "5′ss" == str(scm[5]):
			for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donorGT.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon_inintron.bed"):
				newexon = re.split('[\t\n]',d)
				if str(scm[0]) == str(newexon[0]) and str(scm[1]) == str(newexon[5]) and str(scm[9]) == str(newexon[8]):
					isNothing = False
					outfile.writelines(c[:-1] + "\t" + str(newexon[4]) + "\t" + str(newexon[1]) + "\t" + str(newexon[2]) + "\t" + str(newexon[7]) + "\n")
					break
	if isNothing:
		print (c)
		#outfile.writelines(c[:-1] + "\t" + str(newexon[4]) + "\t" + str(newexon[1]) + "\t" + str(newexon[2]) + "\t" + str(newexon[7]) + "\n")
outfile.close()



num = 0
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_exonlength.gtf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.addexonlength.vcf"):
	scm = re.split('[\t\n]',c)
	exonnum = re.split('[\t\n\s]',str(scm[14]))
	num += 1
	if "In CDS" == str(scm[13]):
		outfile.writelines(">" + str(num) + "\t" + c)
		for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre.gtf"):
			gtf = re.split('[\t\n]',d)
			trans = re.split('[\t\n;]',str(gtf[9]))
			transID =re.split('[\t\n"]',str(trans[1]))
			if str(scm[9]) == str(transID[1]) and str(num) == str(gtf[0]):
				if "CDS" == str(gtf[3]):
					annonum = re.split('[\t\n\s]',str(trans[6]))
					if int(exonnum[2]) > int(annonum[2]):
						exonlength = int(gtf[5]) - int(gtf[4]) + 1
						outfile.writelines(str(gtf[0]) + "\t" + str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" + str(gtf[4]) + "\t" +  str(gtf[5]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + "recaculate" + "\t" + str(exonlength)  + "\t" +str(gtf[9]) + "\n")
					if int(exonnum[2]) == int(annonum[2]):
						outfile.writelines(str(gtf[0]) + "\t" + str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" +  str(scm[-4]) + "\t" +  str(scm[-3]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + "recaculate" + "\t" + str(scm[17]) + "\t" +str(gtf[9]) + "\n")
					if int(exonnum[2]) < int(annonum[2]):
						exonlength = int(gtf[5]) - int(gtf[4]) + 1
						outfile.writelines(str(gtf[0]) + "\t" + str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" + str(gtf[4]) + "\t" +  str(gtf[5]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + "recaculate" + "\t" + str(exonlength)  + "\t" +str(gtf[9]) + "\n")
				if "exon" == str(gtf[3]):
					annonum = re.split('[\t\n\s]',str(trans[6]))
					if int(exonnum[2]) == int(annonum[2]):
						outfile.writelines(str(gtf[0]) + "\t" + str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" + str(scm[-4]) + "\t" +  str(scm[-3]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + str(gtf[8]) + "\t" +str(gtf[9]) + "\n")
					else:
						outfile.writelines(str(gtf[0]) + "\t" + str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" +  str(gtf[4]) + "\t" + str(gtf[5]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + str(gtf[8]) + "\t" + str(gtf[9]) + "\n")
				if "CDS" != str(gtf[3]) and "exon" != str(gtf[3]):	
					 outfile.writelines(str(gtf[0]) + "\t" + str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" +  str(gtf[4]) + "\t" + str(gtf[5]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + str(gtf[8]) + "\t" + str(gtf[9]) + "\n")
					
outfile.close()

#------------------
#change codon gtf
#------------------


	
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_exonlength_pre_countcodon.gtf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_exonlength.gtf"):
        gtf = re.split('[\t\n]',c)
        if ">" in str(gtf[0]):
                amaris = [0]
                atais = [0]
                tameshi ="0"
        else:
                if "CDS" == str(gtf[3]):
                        cds_long = int(gtf[9]) + int(amaris[-1])
                        amari = int(cds_long) % 3
                        if amari == 0:
                                atai = 0
                                amaris.append(amari)
                                atais.append(atai)
                                tameshi += str(atai)
                        if amari == 1:
                                atai = 2
                                amaris.append(amari)
                                atais.append(atai)
                                tameshi += str(atai)
                        if amari == 2:
                                atai = 1
                                amaris.append(amari)
                                atais.append(atai)
                                tameshi += str(atai)
                        outfile.writelines(str(tameshi)  + "\t" + c)
outfile.close()



#-----------------
#extract last CDS
#-----------------
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_exonlength_pre_countcodon_max.gtf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_exonlength_pre_countcodon.gtf"):
	gtf = re.split('[\t\n]',c)
	trans =  re.split('[\t\n;]',gtf[11])
	transID =  re.split('[\t\n"]',trans[1])
	#print (transID)
	gtfnum = re.split('[\t\n"]',trans[6])
	for i in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_lastcds_number.txt"):
		number =  re.split('[\t\n]',i)
		exonnum = re.split('[\t\n\s]',str(number[1]))
		if len(exonnum) == 3:
			if str(number[0]) == str(transID[1]) and str(gtfnum[0]) == str(number[1]):
				outfile.writelines(c)
				break

outfile.close()


#-------------------------
#change codon number
#-------------------------

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_change.gtf","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_exonlength_pre_countcodon_max.gtf"):
	codon = re.split('[\t\n]',i)
	#print (codon)
	count = -1
	for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_exonlength.gtf"):
		gtf = re.split('[\t\n]',c)
		if ">" in str(gtf[0]):
			pass
		else:
			if "CDS" == str(gtf[3]):
				if int(codon[1]) == int(gtf[0]):
					count += 1
					num = codon[0][count]
					outfile.writelines(str(gtf[0]) + "\t" + str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" + str(gtf[4]) + "\t" + str(gtf[5]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + str(num) + "\t" + str(gtf[10]) + "\n")
			if "exon" == str(gtf[3]) or "start_codon" == str(gtf[3]) or "stop_codon"  == str(gtf[3]) or "UTR" == str(gtf[3]) or "transcript" == str(gtf[3]):
					outfile.writelines(str(gtf[0]) + "\t" + str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" + str(gtf[4]) + "\t" + str(gtf[5]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + str(gtf[8]) + "\t" + str(gtf[9]) + "\n")




outfile.close()


#---------------
#uniq
#--------------
cmd = "less -S /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_pre_change.gtf|sort -k1 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_change.gtf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)



try:
	cmd = "mkdir /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)
except:
	pass


#------------------
#Split number
#-----------------
cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_change.gtf|tail -n1|cut -f1 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/max_number"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

'''
split_num = 0
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/max_number"):
	num  = re.split('[\t\n]',c)
	split_num = int(num[0])

split_number = range(split_num + 1)

for d in list(split_number):
	outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(d) + "_scm_change.gtf","w")
	for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/genemap2.Phenotype3_gencode_v29_" + chr + "_scm_change.gtf"):
		gtf = re.split('[\t\n]',c)
		if str(gtf[0]) == str(d):
			outfile.writelines(str(gtf[1]) + "\t" + str(gtf[2]) + "\t" + str(gtf[3]) + "\t" + str(gtf[4]) + "\t" + str(gtf[5]) + "\t" + str(gtf[6]) + "\t" + str(gtf[7]) + "\t" + str(gtf[8]) + "\t" + str(gtf[9]) + "\n")


outfile.close()


#----------
#gffread
#----------


for d in list(split_number):	
	cmd = "gffread -E /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(d) + "_scm_change.gtf -g /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa -o /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(d) + "_scm_change_gffread.gff3"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)	

	cmd = "gffread /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(d) + "_scm_change_gffread.gff3 -g /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/GENCODE/genome/GRCh38/GRCh38.p12.genome.fa -w /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(d) + "_scm_change_gffread.transcripts.fa -x /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(d) + "_scm_change_gffread.cds.fa -y /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gffread/" + str(d) + "_scm_change_gffread.protein.fa"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)
