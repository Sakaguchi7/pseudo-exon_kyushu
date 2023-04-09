#!/usr/bin/python

import sys
import re
import subprocess
import os
import gzip

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
#---------------------------
#Download data & select SNV(ALL finish)
#--------------------------
cmd1 = "python /mnt/data6/narumi/gnomAD/grch38/chr_add_grch38.py"


#------------------------------
#Intersect OMIM gene & SNVs
#-------------------------------
cmd2 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/intersect_omim_gencode_v29.py " + chr 
print (cmd2)
contents = subprocess.check_call(cmd2,shell=True)

#-------------------------
#Select A or G or T base
#-------------------------
cmd3 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/select_AGT_direction_bar_transcript_spt_gencode_v29.py " + chr
print (cmd3)
contents = subprocess.check_call(cmd3,shell=True)

#-------------------------------------------------
#Select acceptor site(chr21 only)
#-------------------------------------------------
cmd4 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/select_acc_region_P3_v2.py " + chr
print (cmd4)
contents = subprocess.check_call(cmd4,shell=True)


#-------------------------------------------------
#Select donor site
#Caculate maxEntScan score of 5'ss (Ref) sequence
#-------------------------------------------------
cmd5 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/select_dono_region_P3_v2.py " +chr
print (cmd5)
contents = subprocess.check_call(cmd5,shell=True)

#----------------------------------------------------
#Caculate maxEntScan score of 3 & 5'ss sequence
#----------------------------------------------------
cmd6 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/change_fasta_snv_before_MES_P3.py " + chr
print (cmd6)
contents = subprocess.check_call(cmd6,shell=True)


#------------------
#MES score ≥ 0
#------------------
cmd7 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/MES_score_over0.py " + chr
print (cmd7)
contents = subprocess.check_call(cmd7,shell=True)



#-------------------------------------------
#make SCM lists
#-------------------------------------------
cmd8 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/temporary_SCMlist_v3.py " + chr
print (cmd8)
contents = subprocess.check_call(cmd8,shell=True)



#-----------------
#Select VCF data
#-----------------

cmd9 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/make_spliceaiformat.py " + chr
print (cmd9)
contents = subprocess.check_call(cmd9,shell=True)


#-----------------
#AF ≤ 0.01
#-----------------
cmd10 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/confirm_AF.py " + chr
print (cmd10)
contents = subprocess.check_call(cmd10,shell=True)



#----------------------
#Separate near exon(50nt)
#----------------------
cmd11 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/separate_region.py " + chr
print (cmd11)
contents = subprocess.check_call(cmd11,shell=True)




#-------------------------
#Near exonic region

#-----------------
#SpliceAI
#-----------------
cmd11 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/spliceAI_v2.py " + chr
print (cmd11)
contents = subprocess.check_call(cmd11,shell=True)



#---------------------
#Δscore ≥ 0.8 or 0.5
#---------------------
cmd12 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/spliceai_score.py " + chr
print (cmd12)
contents = subprocess.check_call(cmd12,shell=True)



#-------------------------------------------
#SNVs in exonic regions or intronic regions
#decide splice site & altexon
#-------------------------------------------
cmd13 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/search_exon_refSS_v3.py " + chr
print (cmd13)
contents = subprocess.check_call(cmd13,shell=True)

#--------------------
#Add ClinVar
#---------------------

cmd14 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/add_clinvar.py " + chr
print (cmd14)
contents = subprocess.check_call(cmd14,shell=True)

#--------------------
#Add ANNOVAR(exon)
#---------------------

cmd15 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/add_annovar.py " + chr
print (cmd15)
contents = subprocess.check_call(cmd15,shell=True)

cmd16 = "python /mnt/data6/narumi/gnomAD/grch38/GTEx/ChrIDwork/GTEx_SCManalysis.py " + chr
print (cmd16)
contents = subprocess.check_call(cmd16,shell=True)




#--------------------------------
#decide splice site & altintron
#----------------------------------

cmd17 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/search_exon_refSS_SNV_intronnumber.py " + chr
print (cmd17)
contents = subprocess.check_call(cmd17,shell=True)

cmd18 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/caculate_length_inintronSNV.py " + chr
print (cmd18)
contents = subprocess.check_call(cmd18,shell=True)


#------------------
#search PTC
#------------------
cmd19 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/search_PTC_v2.py " + chr
print (cmd19)
contents = subprocess.check_call(cmd19,shell=True)

cmd20 = "python /mnt/data6/narumi/gnomAD/grch38/ChrIDwork/find_NMD.py " + chr
print (cmd20)
contents = subprocess.check_call(cmd20,shell=True)


#-------------------
#Caculate PSI value
#-------------------

cmd21 = "python /mnt/data6/narumi/gnomAD/grch38/GTEx/ChrIDwork/select_junction.py " + chr
print (cmd21)
contents = subprocess.check_call(cmd21,shell=True)

'''

#-----------------
#Domain
#----------------

cmd22 = "python /mnt/data6/narumi/OMIM/gffread/domain_protein_sequence_v2.py " + chr
print (cmd22)
contents = subprocess.check_call(cmd22,shell=True)



#refresh_error.py"
