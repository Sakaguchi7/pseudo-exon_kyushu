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

#---------------------
# Acceptor site(A)
#---------------------

print (chr)
print ("accA")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf","w")
outfile3 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed","w")
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_refexon_prefasta.bed","w")
    
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    bed = re.split('[\t\n]',c)
    isNothing = True
    for i in open("/mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + ".gtf"):
        exon =  re.split('[\t\n]',i)
        deta = re.split('[\t;]',str(exon[8]))
        transIDs = re.split('[\t"]',str(deta[1]))
        transcript_ID = str(transIDs[1])
        exon_number = str(deta[6])
        if str(exon[0]) == str(bed[0]) and int(exon[3]) <= int(bed[1]) and int(exon[4]) >= int(bed[1]) and str(bed[9]) == str(transcript_ID):
            isNothing = False
            if "+" == str(exon[6]):
                exon_start = int(bed[1]) + 2
                ref_acc_start = int(exon[3]) - 21
                ref_acc_end = int(exon[3]) + 2
                length = int(exon[4]) - int(exon_start) + 1
                outfile.writelines(c) 
                outfile2.writelines(str(exon[0]) + "\t"  + str(ref_acc_start) + "\t" + str(ref_acc_end) + "\t" + str(transcript_ID) + "\t" + str(exon_number) + "\t" + str(exon[6]) + "\n")
                outfile3.writelines(str(exon[0]) + "\t" + str(exon_start) + "\t" + str(exon[4]) + "\t" + str(exon_number) + "\t" + str(bed[1]) + "\t" + str(exon[6]) + "\t" + str(length) + "\t" + str(transcript_ID) + "\n")
                break
            if "-" == str(exon[6]):
                exon_end = int(bed[1]) - 2
                ref_acc_start = int(exon[4]) - 3
                ref_acc_end = int(exon[4]) + 20
                length = int(exon_end) - int(exon[3]) + 1 
                outfile.writelines(c)
                outfile2.writelines(str(exon[0]) + "\t"  + str(ref_acc_start) + "\t" + str(ref_acc_end) + "\t" + str(transcript_ID) + "\t" + str(exon_number) + "\t" + str(exon[6]) + "\n")
                outfile3.writelines(str(exon[0]) + "\t" + str(exon[3]) + "\t" + str(exon_end) + "\t" + str(exon_number) + "\t" + str(bed[1]) + "\t" + str(exon[6]) + "\t" + str(length) + "\t" + str(transcript_ID) + "\n")
                break



outfile.close()
outfile2.close()
outfile3.close()

outfile4 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    isNothing = True
    for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf"):
        if str(c) == str(i):
            isNothing = False
            break
    if isNothing:
        outfile4.writelines(c)
outfile4.close()

#---------------------
# Acceptor site(G)
#---------------------
print ("accG")

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf","w")
outfile3 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed","w")
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_refexon_prefasta.bed","w")

for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    bed = re.split('[\t\n]',c)
    isNothing = True
    for i in open("/mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + ".gtf"):
        exon =  re.split('[\t\n]',i)
        deta = re.split('[\t;]',str(exon[8]))
        transIDs = re.split('[\t"]',str(deta[1]))
        transcript_ID = str(transIDs[1])
        exon_number = str(deta[6])
        if str(exon[0]) == str(bed[0]) and int(exon[3]) <= int(bed[1]) and int(exon[4]) >= int(bed[1]) and str(bed[9]) == str(transcript_ID):
            isNothing = False
            outfile.writelines(c)
            if "+" == str(exon[6]):
                exon_start = int(bed[1]) + 1
                ref_acc_start = int(exon[3]) - 21
                ref_acc_end = int(exon[3]) + 2
                length = int(exon[4]) - int(exon_start) + 1
                outfile2.writelines(str(exon[0]) + "\t"  + str(ref_acc_start) + "\t" + str(ref_acc_end) + "\t" + str(transcript_ID) + "\t" + str(exon_number) + "\t" + str(exon[6]) + "\n")
                outfile3.writelines(str(exon[0]) + "\t" + str(exon_start) + "\t" + str(exon[4]) + "\t" + str(exon_number) + "\t" + str(bed[1]) + "\t" + str(exon[6]) + "\t" + str(length) + "\t" + str(transcript_ID) + "\n")
                break
            if "-" == str(exon[6]):
                exon_end = int(bed[1]) - 1
                ref_acc_start = int(exon[4]) - 3
                ref_acc_end = int(exon[4]) + 20
                length = int(exon_end) - int(exon[3]) + 1 
                outfile2.writelines(str(exon[0]) + "\t"  + str(ref_acc_start) + "\t" + str(ref_acc_end) + "\t" + str(transcript_ID) + "\t" + str(exon_number) + "\t" + str(exon[6]) + "\n")
                outfile3.writelines(str(exon[0]) + "\t" + str(exon[3]) + "\t" + str(exon_end) + "\t" + str(exon_number) + "\t" + str(bed[1]) + "\t" + str(exon[6]) + "\t" + str(length) + "\t" + str(transcript_ID) + "\n")
                break



outfile.close()
outfile2.close()
outfile3.close()

outfile4 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    isNothing = True
    for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf"):
        if str(c) == str(i):
            isNothing = False
            break
    if isNothing:
        outfile4.writelines(c)
outfile4.close()



#---------------------
# Donor site(G)
#---------------------
print ("donoG")

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf","w")
outfile3 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed","w")
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_refexon_prefasta.bed","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    bed = re.split('[\t\n]',c)
    isNothing = True
    for i in open("/mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + ".gtf"):
        exon =  re.split('[\t\n]',i)
        deta = re.split('[\t;]',str(exon[8]))
        transIDs = re.split('[\t"]',str(deta[1]))
        transcript_ID = str(transIDs[1])
        exon_number = str(deta[6])
        if str(exon[0]) == str(bed[0]) and int(exon[3]) <= int(bed[1]) and int(exon[4]) >= int(bed[1]) and str(bed[9]) == str(transcript_ID):
            isNothing = False
            if "+" == str(exon[6]):
                exon_end = int(bed[1]) - 1
                ref_acc_start = int(exon[4]) - 3
                ref_acc_end = int(exon[4]) + 6
                length = int(exon_end) - int(exon[3]) + 1  
                outfile.writelines(c)
                outfile2.writelines(str(exon[0]) + "\t"  + str(ref_acc_start) + "\t" + str(ref_acc_end) + "\t" + str(transcript_ID) + "\t" + str(exon_number) + "\t" + str(exon[6]) + "\n")
                outfile3.writelines(str(exon[0]) + "\t" + str(exon[3]) + "\t" + str(exon_end) + "\t" + str(exon_number) + "\t" + str(bed[1]) + "\t" + str(exon[6]) + "\t" + str(length) + "\t" + str(transcript_ID) + "\n")
                break
            if "-" == str(exon[6]):
                exon_start = int(bed[1]) + 1
                ref_acc_start = int(exon[3]) - 7
                ref_acc_end = int(exon[3]) + 2
                length = int(exon[4]) - int(exon_start) + 1
                outfile.writelines(c)
                outfile2.writelines(str(exon[0]) + "\t"  + str(ref_acc_start) + "\t" + str(ref_acc_end) + "\t" + str(transcript_ID) + "\t" + str(exon_number) + "\t" + str(exon[6]) + "\n")
                outfile3.writelines(str(exon[0]) + "\t" + str(exon_start) + "\t" + str(exon[4]) + "\t" + str(exon_number) + "\t" + str(bed[1]) + "\t" + str(exon[6]) + "\t" + str(length) + "\t" + str(transcript_ID) + "\n")
                break



outfile.close()
outfile2.close()
outfile3.close()

outfile4 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    isNothing = True
    for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf"):
        if str(c) == str(i):
            isNothing = False
            break
    if isNothing:
        outfile4.writelines(c)
outfile4.close()


#---------------------
# Donor site(T)
#---------------------
print ("donoT")


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf","w")
outfile3 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_newexon.bed","w")
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_SCM_refexon_prefasta.bed","w")

for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    bed = re.split('[\t\n]',c)
    isNothing = True
    for i in open("/mnt/data6/narumi/OMIM/v29_exon/genemap2.Phenotype3_gencode_v29_exon_" + chr + ".gtf"):
        exon =  re.split('[\t\n]',i)
        deta = re.split('[\t;]',str(exon[8]))
        transIDs = re.split('[\t"]',str(deta[1]))
        transcript_ID = str(transIDs[1])
        exon_number = str(deta[6])
        if str(exon[0]) == str(bed[0]) and int(exon[3]) <= int(bed[1]) and int(exon[4]) >= int(bed[1]) and str(bed[9]) == str(transcript_ID):
            isNothing = False
            if "+" == str(exon[6]):
                exon_end = int(bed[1]) - 2
                ref_acc_start = int(exon[4]) - 3
                ref_acc_end = int(exon[4]) + 6
                length = int(exon_end) - int(exon[3]) + 1  
                outfile.writelines(c)
                outfile2.writelines(str(exon[0]) + "\t"  + str(ref_acc_start) + "\t" + str(ref_acc_end) + "\t" + str(transcript_ID) + "\t" + str(exon_number) + "\t" + str(exon[6]) + "\n")
                outfile3.writelines(str(exon[0]) + "\t" + str(exon[3]) + "\t" + str(exon_end) + "\t" + str(exon_number) + "\t" + str(bed[1]) + "\t" + str(exon[6]) + "\t" + str(length) + "\t" + str(transcript_ID) + "\n")
                break
            if "-" == str(exon[6]):
                exon_start = int(bed[1]) + 2
                ref_acc_start = int(exon[3]) - 7
                ref_acc_end = int(exon[3]) + 2
                length = int(exon[4]) - int(exon_start) + 1
                outfile.writelines(c)
                outfile2.writelines(str(exon[0]) + "\t"  + str(ref_acc_start) + "\t" + str(ref_acc_end) + "\t" + str(transcript_ID) + "\t" + str(exon_number) + "\t" + str(exon[6]) + "\n")
                outfile3.writelines(str(exon[0]) + "\t" + str(exon_start) + "\t" + str(exon[4]) + "\t" + str(exon_number) + "\t" + str(bed[1]) + "\t" + str(exon[6]) + "\t" + str(length) + "\t" + str(transcript_ID) + "\n")
                break



outfile.close()
outfile2.close()
outfile3.close()

outfile4 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    isNothing = True
    for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf"):
        if str(c) == str(i):
            isNothing = False
            break
    if isNothing:
        outfile4.writelines(c)
outfile4.close()




