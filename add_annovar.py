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



cmd = "cat /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accA.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_accG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoG.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_donoT.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exon.vcf |cut -f1,2,3,4,5,6,7,8 > /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exoncut.vcf"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)



cmd = "perl /home/narumi/install/annovar/table_annovar.pl /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_exoncut.vcf /home/narumi/install/annovar/humandb_hg38/ -buildver hg38 -out /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/annovar_" + chr + " -remove -protocol dbnsfp35a -operation f -nastring . -vcfinput -polish"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "perl /home/narumi/install/annovar/annotate_variation.pl -out /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/annotate_gencodeV29 -build hg38 /mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/annovar_" + chr + ".avinput /home/narumi/install/annovar/humandb_hg38/ -dbtype wgEncodeGencodeBasicV29"
print (cmd)
contents = subprocess.check_call(cmd,shell=True)


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/annotate_gencodeV29.exonic_variant_function_UTR","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/annovar_" + chr + ".avinput"):
    isNothing = True 
    for d in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/annotate_gencodeV29.exonic_variant_function"):
        if c in d:
            isNothing = False 
            outfile.writelines(d)
            break
    if isNothing:
        outfile.writelines("line0" + "\t" + "UTR" + "\t" + "None" + "\t" + c)

outfile.close()

outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.list.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar.vcf"):
    clinvar = re.split('[\t\n]',c)
    isNothing = True
    for i in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/annotate_gencodeV29.exonic_variant_function_UTR"):
        anno = re.split('[\t\n]',i)
        if str(anno[11]) == str(clinvar[0]) and str(anno[12]) == str(clinvar[1]) and str(anno[13]) == str(clinvar[2]) and str(anno[14]) == str(clinvar[3]) and str(anno[15]) == str(clinvar[4]):
            isNothing = False
            outfile.writelines(c[:-1] + "\t" + str(anno[1]) + "\n")
            break
    if isNothing:
        outfile.writelines(c[:-1] + "\t" + "In intronic" + "\n")

outfile.close()


outfile = open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.listUTR.vcf","w")
for c in open("/mnt/data6/narumi/gnomAD_analysis/grch38/" + chr + "/gnomad.genomes.r3.0.sites." + chr + "_oCDS_P3_allcat.getfasta_v2_all_MESscore_over0_AF0.01h_exonslop50h.spliceai_over08_clinvar_annovar.list.vcf"):
	lists = re.split('[\t\n]',c)
	isNothing = True
	if "UTR" == str(lists[12]):
		isNothing = True
		for i in open("/mnt/data6/narumi/OMIM/genemap2.Phenotype3_gencode_v29_start_codon.gtf"):
			codon = re.split('[\t\n]',i)
			if str(lists[9]) in str(codon[8]):
				isNothing = False
				if "+" == str(codon[6]):
					if int(lists[1]) < int(codon[3]):
						outfile.writelines(str(lists[0]) + "\t" + str(lists[1]) + "\t" + str(lists[2]) + "\t" + str(lists[3]) + "\t" + str(lists[4]) + "\t" + str(lists[5]) + "\t" + str(lists[6]) + "\t" + str(lists[7]) + "\t" + str(lists[8]) + "\t" + str(lists[9]) + "\t" + str(lists[10]) + "\t" + str(lists[11]) + "\t" + "5′UTR" + "\n")
						break
					if int(lists[1]) > int(codon[3]):
						outfile.writelines(str(lists[0]) + "\t" + str(lists[1]) + "\t" + str(lists[2]) + "\t" + str(lists[3]) + "\t" + str(lists[4]) + "\t" + str(lists[5]) + "\t" + str(lists[6]) + "\t" + str(lists[7]) + "\t" + str(lists[8]) + "\t" + str(lists[9]) + "\t" + str(lists[10]) + "\t" + str(lists[11]) + "\t" + "3′UTR" + "\n")
						break
				if "-" == str(codon[6]):
					if int(lists[1]) < int(codon[3]):
						outfile.writelines(str(lists[0]) + "\t" + str(lists[1]) + "\t" + str(lists[2]) + "\t" + str(lists[3]) + "\t" + str(lists[4]) + "\t" + str(lists[5]) + "\t" + str(lists[6]) + "\t" + str(lists[7]) + "\t" + str(lists[8]) + "\t" + str(lists[9]) + "\t" + str(lists[10]) + "\t" + str(lists[11]) + "\t" + "3′UTR" + "\n")
						break
					if int(lists[1]) > int(codon[3]):
						outfile.writelines(str(lists[0]) + "\t" + str(lists[1]) + "\t" + str(lists[2]) + "\t" + str(lists[3]) + "\t" + str(lists[4]) + "\t" + str(lists[5]) + "\t" + str(lists[6]) + "\t" + str(lists[7]) + "\t" + str(lists[8]) + "\t" + str(lists[9]) + "\t" + str(lists[10]) + "\t" + str(lists[11]) + "\t" + "5′UTR" + "\n")
						break
		if isNothing:
			outfile.writelines(str(lists[0]) + "\t" + str(lists[1]) + "\t" + str(lists[2]) + "\t" + str(lists[3]) + "\t" + str(lists[4]) + "\t" + str(lists[5]) + "\t" + str(lists[6]) + "\t" + str(lists[7]) + "\t" + str(lists[8]) + "\t" + str(lists[9]) + "\t" + str(lists[10]) + "\t" + str(lists[11]) + "\t" + "lincRNA" + "\n")
	else:
		isNothing = True
		outfile.writelines(c)

outfile.close()
