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



#--------------------
#Acceptor site (G)
#--------------------

print ("accG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2.bed"):
        snp = re.split('[\t\n]',i)
        ennki = list(str(snp[6]))
        if "+" == str(snp[5]) and "A" == str(ennki[18]):
                if str(snp[3]) == str(ennki[19]) :
                        ennki[19] = str(snp[4])
                       	after = str(''.join(ennki))
                       	outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
        if "-" == str(snp[5]) and "A" == str(ennki[18]):
                if str(snp[3]) == "A" and "T" == str(ennki[19]) :
                        ennki[19] = "G"
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
                if str(snp[3]) == "T" and "A" == str(ennki[19]) :
                        ennki[19] = "G"
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
                if str(snp[3]) == "G" and "C" == str(ennki[19]) :
                        ennki[19] = "G"
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
	
outfile.close()
#----------------------
#Acceptor site (A)
#------------------------

print ("accA")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2.bed"):
        snp = re.split('[\t\n]',i)
        ennki = list(str(snp[6]))
        if "+" == str(snp[5]) and "G" == str(ennki[19]):
                if str(snp[3]) == str(ennki[18]):
                        ennki[18] = str(snp[4])
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
        if "-" == str(snp[5]) and "G" == str(ennki[19]):
                if str(snp[3]) == "A" and "T" == str(ennki[18]) :
                        ennki[18] = "A"
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
                if str(snp[3]) == "G" and "C" == str(ennki[18]) :
                        ennki[18] = "A"
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
                if str(snp[3]) == "C" and "G" == str(ennki[18]) :
                        ennki[18] = "A"
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
outfile.close()

print ("donoT")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2.bed"):
        snp = re.split('[\t\n]',i)
        ennki = list(str(snp[6]))
        if "+" == str(snp[5]) and "G" == str(ennki[3]):
                if str(snp[3]) == str(ennki[4]):
                        ennki[4] = str(snp[4])
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
        if "-" == str(snp[5]) and "G" == str(ennki[3]):
                if str(snp[3]) == "T" and "A" == str(ennki[4]) :
                        ennki[4] = "T" 
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
                if str(snp[3]) == "G" and "C" == str(ennki[4]) :
                        ennki[4] = "T" 
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
                if str(snp[3]) == "C" and "G" == str(ennki[4]) :
                        ennki[4] = "T" 
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
        
outfile.close()

print ("donoG")
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt.bed","w")
for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2.bed"):
        snp = re.split('[\t\n]',i)
        ennki = list(str(snp[6]))
        if "+" == str(snp[5]) and "T" == str(ennki[4]):
                if str(snp[3]) == str(ennki[3]):
                        ennki[3] = str(snp[4])
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
        if "-" == str(snp[5]) and "T" == str(ennki[4]): 
                if str(snp[3]) == "G" and "C" == str(ennki[3]) :
                        ennki[3] = "G" 
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
                if str(snp[3]) == "T" and "A" == str(ennki[3]) :
                        ennki[3] = "G" 
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
                if str(snp[3]) == "A" and "T" == str(ennki[3]) :
                        ennki[3] = "G" 
                        after = str(''.join(ennki))
                        outfile.writelines(str(i[:-1]) + "\t" + str(after) + "\n")
outfile.close()

#------------------------
#Caculate MES Score
#------------------------

#----------------------
#Acceptor site (A)
#------------------------
cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt.bed|cut -f1,2,3,4,5,6 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt_cut.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt.bed|cut -f 7 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_ref.fasta"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt.bed|cut -f 8 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt.fasta"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "perl score3.pl /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_ref.fasta > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_ref.maxentscan.score"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "perl score3.pl /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt.fasta > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt.maxentscan.score"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "paste /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt_cut.bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_ref.maxentscan.score /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_alt.maxentscan.score > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

#----------------------
#Acceptor site (G)
#------------------------

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt.bed|cut -f1,2,3,4,5,6 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt_cut.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt.bed|cut -f 7 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_ref.fasta"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt.bed|cut -f 8 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt.fasta"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "perl score3.pl /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_ref.fasta > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_ref.maxentscan.score"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "perl score3.pl /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt.fasta > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt.maxentscan.score"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)


cmd = "paste /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt_cut.bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_ref.maxentscan.score /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_alt.maxentscan.score > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)
	


#----------------------
#Donor site (G)
#------------------------

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt.bed|cut -f1,2,3,4,5,6 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt_cut.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt.bed|cut -f 7 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_ref.fasta"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt.bed|cut -f 8 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt.fasta"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "perl score5.pl /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_ref.fasta > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_ref.maxentscan.score"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "perl score5.pl /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt.fasta > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt.maxentscan.score"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)



cmd = "paste /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt_cut.bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_ref.maxentscan.score /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_alt.maxentscan.score > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)
#----------------------
#Donor site (T)
#------------------------

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt.bed|cut -f1,2,3,4,5,6 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt_cut.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt.bed|cut -f 7 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_ref.fasta"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "less /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt.bed|cut -f 8 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt.fasta"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "perl score5.pl /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_ref.fasta > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_ref.maxentscan.score"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "perl score5.pl /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt.fasta > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt.maxentscan.score"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "paste /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt_cut.bed /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_ref.maxentscan.score /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_alt.maxentscan.score > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore.bed"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

