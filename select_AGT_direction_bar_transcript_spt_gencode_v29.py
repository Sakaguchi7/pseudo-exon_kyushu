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

print (chr)
outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.vcf","w")
outfile2 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.vcf","w")
outfile3 = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoaccG.vcf","w")
f =  open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_snv_omimCDS_Phenotype3.vcf")
pbar = tqdm(total=len(open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr+ "/gnomad.genomes.r3.0.sites." + chr + "_snv_omimCDS_Phenotype3.vcf").readlines()))
for i in tqdm(f):
	pbar.update(1)
	if "#" in i:
		outfile.writelines(i)
		outfile2.writelines(i)
		outfile3.writelines(i)
	else:
		snp = re.split('[\t\n]',i)
		if len(snp[3]) == 1 and len(snp[4]) == 1:
			for c in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_gene.gtf"):
				gtf = re.split('[\t\n]',c)
				if str(snp[0]) == str(gtf[0]) and int(snp[1]) >= int(gtf[3]) and int(snp[1]) <= int(gtf[4]):
					if "+" == str(gtf[6]):
						if "A" == str(snp[4]):
							outfile.writelines(i[:-1] + "\t" + str(gtf[6]) + "\n")
							break
						if "T" == str(snp[4]):
							outfile2.writelines(i[:-1] + "\t" + str(gtf[6]) + "\n")
							break
						if "G" == str(snp[4]):
							outfile3.writelines(i[:-1] + "\t" + str(gtf[6]) + "\n")
							break
					if "-" == str(gtf[6]):
						if "T" == str(snp[4]):
							outfile.writelines(i[:-1] + "\t" + str(gtf[6]) + "\n")
							break
						if "A" == str(snp[4]):
							outfile2.writelines(i[:-1] + "\t" + str(gtf[6]) + "\n")
							break
						if "C" == str(snp[4]):
							outfile3.writelines(i[:-1] + "\t" + str(gtf[6]) + "\n")
							break
							
outfile.close()
outfile2.close()
outfile3.close()
pbar.close()
