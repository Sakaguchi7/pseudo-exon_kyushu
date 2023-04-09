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
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    scms = re.split('[\t\n]',c)
    ACS = re.split('[\t\n;]',str(scms[7]))
    isNothing = True
    for i in open("/mnt/data6/narumi/clinvar/clinvar_20200919_" + chr + ".vcf"):
        if "#" in i:
            pass
        else:
            clinvar = re.split('[\t\n]',i)
            #detail = re.split('[\t\n]',str(clinvar[7])) 
        #print (clinvar)
            if "chr" + str(clinvar[0]) == str(scms[0]) and str(clinvar[1]) == str(scms[1]) and str(clinvar[3]) == str(scms[3]) and str(clinvar[4]) == str(scms[4]) :
                isNothing = False
                #print (detail[)
                numbers = str(clinvar[7]).find('CLNSIG=')
                if len(str(numbers)) > 2:
                #print (str(clinvar[7]))
                        det = str(clinvar[7])[int(numbers):-1]
               	        detail = re.split('[\t\n;=]',str(det))
                #variants = re.split('[\t\n|]',str(det))
                #print (clinvar[7])
                #print (variants)
                        outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "3′ss" + "\t" + str(ACS[0]) + "\t" + str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + str(detail[1]) + "\n")
                        break
                if len(str(numbers)) ==  1:
                        outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "3′ss" + "\t" + str(ACS[0]) + "\t" +     str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + "Check ClinVar" + "\n")
    if isNothing:
        outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "3′ss" + "\t" + str(ACS[0]) + "\t" + str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + "Not in ClinVar" + "\n")

outfile.close()


#---------------------
# Acceptor site(G)
#---------------------
print ("accG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    scms = re.split('[\t\n]',c)
    ACS = re.split('[\t\n;]',str(scms[7]))
    isNothing = True
    for i in open("/mnt/data6/narumi/clinvar/clinvar_20200919_" + chr + ".vcf"):
        if "#" in i:
            pass
        else:
            clinvar = re.split('[\t\n]',i)
            #detail = re.split('[\t\n]',str(clinvar[7])) 
        #print (clinvar)
            if "chr" + str(clinvar[0]) == str(scms[0]) and str(clinvar[1]) == str(scms[1]) and str(clinvar[3]) == str(scms[3]) and str(clinvar[4]) == str(scms[4]) :
                isNothing = False
                #print (detail[)
                numbers = str(clinvar[7]).find('CLNSIG=')
                if len(str(numbers)) > 2:
                #print (str(clinvar[7]))
                         det = str(clinvar[7])[int(numbers):-1]
                         detail = re.split('[\t\n;=]',str(det))
                         outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "3′ss" + "\t" + str(ACS[0]) + "\t" + str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + str(detail[1]) + "\n")
                         break
                if len(str(numbers)) ==  1:
                         outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "3′ss" + "\t" + str(ACS[0]) + "\t" +     str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + "Check ClinVar" + "\n")
    if isNothing:
        outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "3′ss" + "\t" + str(ACS[0]) + "\t" + str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + "Not in ClinVar" + "\n")

outfile.close()

#---------------------
# Donor site(G)
#---------------------
print ("donoG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    scms = re.split('[\t\n]',c)
    ACS = re.split('[\t\n;]',str(scms[7]))
    isNothing = True
    for i in open("/mnt/data6/narumi/clinvar/clinvar_20200919_" + chr + ".vcf"):
        if "#" in i:
            pass
        else:
            clinvar = re.split('[\t\n]',i)
            #detail = re.split('[\t\n]',str(clinvar[7])) 
        #print (clinvar)
            if "chr" + str(clinvar[0]) == str(scms[0]) and str(clinvar[1]) == str(scms[1]) and str(clinvar[3]) == str(scms[3]) and str(clinvar[4]) == str(scms[4]) :
                isNothing = False
                #print (detail[)
                numbers = str(clinvar[7]).find('CLNSIG=')
                if len(str(numbers)) > 2:
                #print (str(clinvar[7]))
                        det = str(clinvar[7])[int(numbers):-1]
                       	#print (det)
                        detail = re.split('[\t\n;=]',str(det))
                        #print (scms)
                        #print (ACS)
                #print (det)
                        outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "5′ss" + "\t" + str(ACS[0]) + "\t" + str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + str(detail[1]) + "\n")
                        break
                if len(str(numbers)) ==  1:
                        outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "5′ss" + "\t" + str(ACS[0]) + "\t" + str(ACS[1]) + "\t" +str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" +  "Check ClinVar" + "\n")
    if isNothing:
        outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "5′ss" + "\t" + str(ACS[0]) + "\t" + str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + "Not in ClinVar" + "\n")

outfile.close()

#---------------------
# Donor site(T)
#---------------------
print ("donoT")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08.vcf"):
    scms = re.split('[\t\n]',c)
    ACS = re.split('[\t\n;]',str(scms[7]))
    isNothing = True
    for i in open("/mnt/data6/narumi/clinvar/clinvar_20200919_" + chr + ".vcf"):
        if "#" in i:
            pass
        else:
            clinvar = re.split('[\t\n]',i)
            #detail = re.split('[\t\n]',str(clinvar[7])) 
        #print (clinvar)
            if "chr" + str(clinvar[0]) == str(scms[0]) and str(clinvar[1]) == str(scms[1]) and str(clinvar[3]) == str(scms[3]) and str(clinvar[4]) == str(scms[4]) :
                isNothing = False
                #print (detail[)
                numbers = str(clinvar[7]).find('CLNSIG=')
                #print (str(clinvar[7]))
                det = str(clinvar[7])[int(numbers):-1]
                #print (det)
                detail = re.split('[\t\n;=]',str(det))
                outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "5′ss" + "\t" + str(ACS[0]) + "\t" + str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + str(detail[1]) + "\n")
                break
    if isNothing:
        outfile.writelines(str(scms[0]) + "\t" + str(scms[1]) + "\t" + str(scms[2]) + "\t" + str(scms[3]) + "\t" + str(scms[4]) + "\t" + "5′ss" + "\t" + str(ACS[0]) + "\t" + str(ACS[1]) + "\t" + str(ACS[2]) + "\t" + str(scms[9]) + "\t" + str(scms[8]) + "\t" + "Not in ClinVar" + "\n")

outfile.close()

cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf |sort -k1,2 -V  > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

