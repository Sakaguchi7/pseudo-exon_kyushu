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
print ("accA")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron_exonnumber.vcf","w")

for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron.vcf"):
    bed = re.split('[\t\n]',c)
    isNothing = True 
    for i in open("/mnt/data6/narumi/OMIM/v29_intron/genemap2.Phenotype3_gencode_v29_intronnumber_" + chr +  ".gtf"):
        exon =  re.split('[\t\n|]',i)
        num = re.split('[\t\s]',str(exon[7]))
        transcript_ID = str(exon[4]) 
        if str(exon[0]) == str(bed[0]) and int(exon[1]) < int(bed[1]) and int(exon[2]) >= int(bed[1]) and str(bed[9]) == str(transcript_ID):
            isNothing = False
            if "+" == str(exon[6]):
                minum = int(num[1]) + 1
                outfile.writelines(c[:-1] + "\t" + str(exon[7]) + "\t" + "exon_number " + str(minum) + "\n")
                break
            if "-" == str(exon[6]):
                minum = int(num[1]) + 1
                outfile.writelines(c[:-1] + "\t" + str(exon[7]) + "\t" + "exon_number " + str(minum) + "\n") 
                break
    if isNothing:
        outfile.writelines(c[:-1] + "\t" + "intron_number unknown" + "\t" + "exon_number unknown"+ "\n")  
outfile.close()


#---------------------
# Acceptor site(G)
#---------------------

print ("accG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron_exonnumber.vcf","w")

for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron.vcf"):
    bed = re.split('[\t\n]',c)
    isNothing = True 
    for i in open("/mnt/data6/narumi/OMIM/v29_intron/genemap2.Phenotype3_gencode_v29_intronnumber_" + chr +  ".gtf"):
        exon =  re.split('[\t\n|]',i)
        num = re.split('[\t\s]',str(exon[7]))
        transcript_ID = str(exon[4]) 
        if str(exon[0]) == str(bed[0]) and int(exon[1]) < int(bed[1]) and int(exon[2]) >= int(bed[1]) and str(bed[9]) == str(transcript_ID):
            isNothing = False
            if "+" == str(exon[6]):
                minum = int(num[1]) + 1
                outfile.writelines(c[:-1] + "\t" + str(exon[7]) + "\t" + "exon_number " + str(minum) + "\n")
                break
            if "-" == str(exon[6]):
                minum = int(num[1]) + 1
                outfile.writelines(c[:-1] + "\t" + str(exon[7]) + "\t" + "exon_number " + str(minum) + "\n") 
                break
    if isNothing:
        outfile.writelines(c[:-1] + "\t" + "intron_number unknown" + "\t" + "exon_number unknown"+ "\n")  
outfile.close()



#---------------------
# Donor site(G)
#---------------------
print ("donoG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron_exonnumber.vcf","w")

for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron.vcf"):
    bed = re.split('[\t\n]',c)
    isNothing = True 
    for i in open("/mnt/data6/narumi/OMIM/v29_intron/genemap2.Phenotype3_gencode_v29_intronnumber_" + chr +  ".gtf"):
        exon =  re.split('[\t\n|]',i)
        num = re.split('[\t\s]',str(exon[7]))
        transcript_ID = str(exon[4]) 
        if str(exon[0]) == str(bed[0]) and int(exon[1]) < int(bed[1]) and int(exon[2]) >= int(bed[1]) and str(bed[9]) == str(transcript_ID):
            isNothing = False
            if "+" == str(exon[6]):
                minum = int(num[1])
                outfile.writelines(c[:-1] + "\t" + str(exon[7]) + "\t" + "exon_number " + str(minum) + "\n")
                break
            if "-" == str(exon[6]):
                minum = int(num[1])
                outfile.writelines(c[:-1] + "\t" + str(exon[7]) + "\t" + "exon_number " + str(minum) + "\n") 
                break
    if isNothing:
        outfile.writelines(c[:-1] + "\t" + "intron_number unknown" + "\t" + "exon_number unknown"+ "\n")  
outfile.close()


#---------------------
# Donor site(T)
#---------------------
print ("donoT")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron_exonnumber.vcf","w")

for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_intron.vcf"):
    bed = re.split('[\t\n]',c)
    isNothing = True 
    for i in open("/mnt/data6/narumi/OMIM/v29_intron/genemap2.Phenotype3_gencode_v29_intronnumber_" + chr +  ".gtf"):
        exon =  re.split('[\t\n|]',i)
        num = re.split('[\t\s]',str(exon[7]))
        transcript_ID = str(exon[4]) 
        if str(exon[0]) == str(bed[0]) and int(exon[1]) < int(bed[1]) and int(exon[2]) >= int(bed[1]) and str(bed[9]) == str(transcript_ID):
            isNothing = False
            if "+" == str(exon[6]):
                minum = int(num[1])
                outfile.writelines(c[:-1] + "\t" + str(exon[7]) + "\t" + "exon_number " + str(minum) + "\n")
                break
            if "-" == str(exon[6]):
                minum = int(num[1])
                outfile.writelines(c[:-1] + "\t" + str(exon[7]) + "\t" + "exon_number " + str(minum) + "\n") 
                break
    if isNothing:
        outfile.writelines(c[:-1] + "\t" + "intron_number unknown" + "\t" + "exon_number unknown"+ "\n")  
outfile.close()

