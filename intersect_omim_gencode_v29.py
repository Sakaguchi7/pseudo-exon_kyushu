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
#------------------------
#contain SNV in OMIM genebody
#------------------------

cmd = "bedtools intersect -a /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_snv.vcf -b /mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_gene.gtf -wa|sort -k1,2 -V |uniq > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_snv_omimCDS_Phenotype3.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)
	


