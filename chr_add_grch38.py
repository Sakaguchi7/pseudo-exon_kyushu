#!/usr/bin/python

import sys
import re
import subprocess
import os






#chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","X"]
#chrs = ["X"]
#chrs = ["22"]
#chrs = ["Y"]

#--------------
#vcf download
#--------------
for chr in chrs:
	chr = "chr" + chr
	cmd = "wget https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites." + str(chr) + ".vcf.bgz"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True) 





#------------------
#change "bgz" to "gz"
#------------------
	cmd = "mv /mnt/data6/narumi/gnomAD/grch38/gnomad.genomes.r3.0.sites." + str(chr) + ".vcf.bgz /mnt/data6/narumi/gnomAD/grch38/gnomad.genomes.r3.0.sites." + str(chr) + ".vcf.gz"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)

	cmd = "gunzip /mnt/data6/narumi/gnomAD/grch38/gnomad.genomes.r3.0.sites." + str(chr) + ".vcf.gz"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)
	
	try:
		cmd = "mkdir /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr
		contents = subprocess.check_call(cmd,shell=True)
	except:
		pass

	outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + str(chr) + "_snv.vcf","w")
	outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + ".header","w")
	for i in open("/mnt/data6/narumi/gnomAD/grch38/gnomad.genomes.r3.0.sites." + str(chr) + ".vcf"):
		if "#" in i:
			outfile.writelines(i)
			outfile2.writelines(i)
		else:
			snp = re.split('[\t\n]',i)
			if len(snp[3]) == 1 and len(snp[4]) == 1:
				outfile.writelines(i)

	outfile.close()
	outfile2.close()


chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
for chr in chrs:
	chr = "chr" + chr
	cmd = "bgzip /mnt/data6/narumi/gnomAD/grch38/gnomad.genomes.r3.0.sites." + str(chr) + ".vcf"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)






'''

#------------------------
#contain SNV in genebody
#------------------------

#	cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + ".vcf -b /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/Ensembl/hg19/ensembl_hg19_start_stop.gtf -wa|sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + "_CDS.vcf"
#	print (cmd)
#	contents = subprocess.check_call(cmd,shell=True)

#-----------------------------
#remove SNV in exonic region
#-----------------------------

#	cmd = "cat /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + ".header /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + "_CDS.vcf > /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + "_CDSh.vcf"
#	print (cmd)
#	contents = subprocess.check_call(cmd,shell=True)

#	cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + "_CDSh.vcf -b /mnt/houman/narumi/gEUVADIS_RNAdata_twins1/Ensembl/hg19/ensembl_hg19_exon.gtf -v > /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites."+ chr + "_CDS_exc_exon.vcf"
#	print (cmd)
#	contents = subprocess.check_call(cmd,shell=True)

	cmd = "bgzip /mnt/data6/narumi/gnomAD/gnomad.genomes.r2.1.1.sites." + num + ".vcf"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)

	cmd = "bgzip /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + ".vcf"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)

	cmd = "bgzip /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + "_CDS.vcf"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)

	cmd = "bgzip /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + "_CDSh.vcf"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)

	
	cmd = "cat /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites." + chr + ".header /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites."+ chr + "_CDS_exc_exon.vcf > /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites."+ chr + "_CDS_exc_exonh.vcf"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)

	cmd = "bgzip /mnt/data6/narumi/gnomAD_analysis/" + chr + "/gnomad.genomes.r2.1.1.sites."+ chr + "_CDS_exc_exon.vcf"
	print (cmd)
	contents = subprocess.check_call(cmd,shell=True)

'''
