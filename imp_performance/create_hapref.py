import re
import sys

#Get individual ID and chromosome
ind = sys.argv[1]
chr = sys.argv[2]

#Define directories and files
out_dir = "/well/hill/adriank/all_imputations/refs_and_vals/references"
hap_dir = "/well/hill/adriank/MEX_BImp/data/1KGP/haps"
hap_file = hap_dir + "/1000GP_Phase3_chr" + chr + ".hap"
amr_ids = "/well/hill/adriank/MEX_BImp/data/array_positions/AMR_ids.txt"
sample = "/well/hill/adriank/MEX_BImp/data/1KGP/1000GP_Phase3.sample"

#Find sample column number
ind_num=0
with open(sample, "r") as samp:
    next(samp)
    for line in samp:
        if not re.search(ind, line):
            ind_num+=1
        else:
            ind_num+=1
            break

#In this case, the haplotye data is two columns for each individual. Since python is base 0 then the first individual
#would be columns 0 and 1. That is 2N -1 and 2N-2.

firsthap = 2*ind_num - 1
secondhap = 2*ind_num - 2

#Find out number of individual in regard of AMR individuals
out_num = 0
with open(amr_ids, "r") as amr:
    for line in amr:
        if not re.search(ind, line):
            out_num+=1
        else:
            out_num+=1
            break
#Define outputs

outhap = out_dir + "/chr" + chr + "/ind_" + out_num + "ref.hap"
outsample = out_dir + "/chr" + chr + "/ind_" + out_num + ".sample"

#Write new haps

#Read hap file
out=open(outhap, "w")
with open(hap_file) as haps:
    for line in haps:
        fields=re.split(" ", line) #split line in individual haps
        del fields[firsthap] #delete first hap the bigger number to not upset the index numbering for the lesser number
        del fields[secondhap] #delete second hap
        out.write(" ".join([str(s) for s in fields]))
out.close()

#Write new sample

outs=open(outsample, "w")
with open(sample, "r") as samp:
    for line in samp:
        if not re.search(ind, line):
            outs.write(line)
        else:
            next
