#!/home/narumi/anaconda3/bin/python
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



cmd = "/mnt/data6/narumi/gnomAD/grch38/GTEx/ChrIDwork/compare_SNVandGTEX_nearexon.py " + chr
print (cmd)
contents = subprocess.check_call(cmd,shell=True)

cmd = "/mnt/data6/narumi/gnomAD/grch38/GTEx/ChrIDwork/find_GTExID_blood_nearexon.py " + chr
print (cmd)
contents = subprocess.check_call(cmd,shell=True)
