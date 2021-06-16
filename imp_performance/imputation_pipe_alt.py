import os
import sys
import time
import configparser
import re
import argparse
import shlex
import subprocess
import pandas as pd

################################################################################
################## READ CONFIG FILE ############################################
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config",dest="config", required=True,help="Path to configuration file")

args = parser.parse_args()
config_file=args.config

################################################################################
############## DEFINE SCRIPT FUNCTIONS #########################################
################################################################################
def al_done(files):
    test = [os.path.exists(f) for f in files]
    if all(test):
        return(True)
    else:
        return(False)

def poll_for_done(flagfile, sleep_time):
    while not os.path.exists(flagfile):
        time.sleep(int(sleep_time))

def poll_for_done_mult(flags, sleep_time):
    finish = False
    while not finish:
        files = [os.path.exists(f) for f in flags]
        if all(files):
            finish = True
        else:
            time.sleep(int(sleep_time))


def write_qsub_script(run_file, job_name, job_project, job_queue, outfile, errfile, run_cmd, flagfile):
    with open(run_file, "w") as script_file:
        script_file.write(
'''#!/bin/bash
#$ -S /bin/bash
#$ -N {0}
#$ -P {1}
#$ -q {2}
#$ -o {3}
#$ -e {4}

{5}

touch {6}
exit 0'''.format(job_name, job_project, job_queue, outfile, errfile, run_cmd, flagfile))

def write_qsub_array_script(run_file,job_name, job_project, job_queue, outfile, errfile, array_start, array_end, run_cmd, flagfile):
    with open(run_file, "w") as script_file:
        script_file.write(
'''#!/bin/bash
#$ -S /bin/bash
#$ -N {0}
#$ -P {1}
#$ -q {2}
#$ -o {3}
#$ -e {4}
#$ -t {5}-{6}

{7}

touch {8}.${{SGE_TASK_ID}}
exit 0'''.format(job_name, job_project, job_queue, outfile, errfile, array_start, array_end, run_cmd, flagfile))


def unphase_plink(plink, file_name):
    command="{0} --vcf {1}_sorted.vcf.gz --keep-allele-order --vcf-half-call m --vcf-filter --make-bed --out {1}_unphased".format(plink, file_name)
    return(command)

def impute_oneref(impute,targ, gen_map, ref_hap, ref_leg,out):
    bounds = "$strt $end"
    command="{0} -filt_rules_l 'TYPE!=Biallelic_SNP' -use_prephased_g -known_haps_g {1} -m {2} -h {3} -l {4} -Ne 20000 -k_hap 2500 -int {5} -o {6}".format(impute, targ, gen_map, ref_hap, ref_leg, bounds, out)
    return(command)

def impute_tworef(impute, targ, gen_map, ref_hap, ref_leg, nat_hap, nat_leg, out):
    bounds = "$strt $end"
    hap_ref = "{0} {1}".format(ref_hap, nat_hap)
    leg_ref = "{0} {1}".format(ref_leg, nat_leg)
    command="{0} -merge_ref_panels -filt_rules_l 'TYPE!=Biallelic_SNP' -use_prephased_g -known_haps_g {1} -m {2} -h {3} -l {4} -Ne 20000 -k_hap 2500 190 -int {5} -o {6}".format(impute, targ, gen_map, hap_ref, leg_ref, bounds, out)
    return(command)

################################################################################
################## GET CONFIG INFO #############################################
################################################################################

#Read config file and get the relevant variables
config = configparser.ConfigParser(allow_no_value=True)
config.read(config_file)

# Get global variables
queue = config.get("Global", "queue name")
project = config.get("Global", "project name")
array_start = int(config.get("Global", "array_start"))
array_end = int(config.get("Global", "array_end"))
log_dir = config.get("Global", "log_dir")
start_chr = int(config.get("Global", "start_chr"))
end_chr = int(config.get("Global", "end_chr"))
step = config.get("Global", "step")
imp_mode = config.get("Global", "mode")

# Get progam paths

vcftools = config.get("Programs", "vcftools")
vcf_sort = config.get("Programs", "vcf-sort")
bcftools= config.get("Programs", "bcftools")
bgzip = config.get("Programs", "bgzip")
tabix = config.get("Programs", "tabix")
perl = config.get("Programs", "perl")
ref_script = config.get("Programs", "ref_script")
val_script = config.get("Programs", "val_script")
shapeit = config.get("Programs", "shapeit")
impute2 = config.get("Programs", "impute2")
R = config.get("Programs", "R")
plink = config.get("Programs", "plink")
bind_script = config.get("Programs", "bind_script")
sums_script = config.get("Programs", "sums_script")
infomaf_script = config.get("Programs", "infomaf_script")
plotmafinfo_script = config.get("Programs", "plotmafinfo_script")
concor_script = config.get("Programs", "concor_script")

# Get files paths

id_file = config.get("Files", "id file")
array = config.get("Files", "array")
vcf_dir = config.get("Files", "vcf_dir")
vcf_name = config.get("Files", "vcf_name")
array_dir = config.get("Files", "array_dir")
map_dir = config.get("Files", "map_dir")
hap_dir = config.get("Files", "hap_dir")
index_dir = config.get("Files", "index_dir")
second_refdir = config.get("Files", "second_refdir")
second_basename = config.get("Files", "second_basename")

#Get output dirs

script_dir = config.get("OutDirs", "script_dir")
out_refdir = config.get("OutDirs", "out_refdir")
val_dir = config.get("OutDirs", "val_dir")
tar_dir = config.get("OutDirs", "tar_dir")
imp_res = config.get("OutDirs", "imp_res")
flag_dir = config.get("OutDirs", "flag_dir")
sums_dir = config.get("OutDirs", "sums_dir")
infomaf_dir = config.get("OutDirs", "infomaf_dir")
figure_dir = config.get("OutDirs", "figure_dir")

#Script names

make_refname = config.get("ScriptNames", "make_refname")
val_tarname = config.get("ScriptNames", "val_tarname")
imp_name = config.get("ScriptNames", "imp_name")
bind_name = config.get("ScriptNames", "bind_name")
sums_name = config.get("ScriptNames", "sums_name")
infomaf_name = config.get("ScriptNames", "infomaf_name")
plotmafinfo_name = config.get("ScriptNames", "plotmafinfo_name")
concor_name = config.get("ScriptNames", "concor_name")

for chrom in range(start_chr,end_chr +1): #Rember that range is not inclusive

#Define some common variables
    vcf_filename = vcf_name.replace("FF", str(chrom))
    vcf_file=os.path.join(vcf_dir,vcf_filename)
    get_id= "ids_file={0}".format(id_file) + "\n" + "id=$(awk \"NR==$SGE_TASK_ID\" $ids_file)"
    gen_map=os.path.join(map_dir, "genetic_map_chr{0}_combined_b37.txt".format(chrom))
    #Define reference legend, hap and sample files for phasing and imputation
    ref_hap=os.path.join(hap_dir, "1000GP_Phase3_chr{0}.hap.gz".format(chrom))
    ref_leg=os.path.join(hap_dir, "1000GP_Phase3_chr{0}.legend.gz".format(chrom))
    ref_sample=os.path.join(hap_dir, "1000GP_Phase3.sample")
    ###########################################################################
    ## CREATE REFERENCE PANEL #################################################
    ###########################################################################

    #Define and create output dir
    out_dir=os.path.join(out_refdir, "chr{0}".format(chrom))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    #Define and create log_dir for reference making
    ref_log=os.path.join(log_dir, "ref_log")
    if not os.path.exists(ref_log):
        os.makedirs(ref_log)

    #Define the base name for the outputs file
    base_name=os.path.join(out_dir, "ind_${SGE_TASK_ID}_ref")


    #Remove individual from the reference
    subset="{0} --gzvcf {1} --remove-indv ${{id}} --recode --recode-INFO-all --stdout | {2} -c > {3}.vcf.gz".format(vcftools, vcf_file, bgzip, base_name)

    #Index sorted vcf file
    index="{0} -p vcf {1}.vcf.gz".format(tabix, base_name)

    #Turn into ref panel
    makeref="{0} {1} -vcf {2}.vcf.gz -leghap {2} -chr {3} -snps_only".format(perl, ref_script, base_name, chrom)

    #Remove unsorted and sorted vcf.
    remove= "rm {0}.vcf.gz*".format(base_name)

    #Join all the commands into a single string
    run_cmd = "\n".join([get_id, subset, index, makeref, remove])

    #Define the script name and write the script
    script_name=os.path.join(script_dir,"{0}_{1}.sh".format(make_refname, chrom))
    flag_file = os.path.join(flag_dir, "chr{0}_flagfile".format(chrom))
    #Check if the script and the out files already exist
    if not os.path.exists(script_name):
        job_name = "{0}_{1}".format(make_refname, chrom)
        logfile = os.path.join(ref_log, "chr{0}_reflog.txt".format(chrom))
        errfile = os.path.join(log_dir, "chr{0}_referr.txt".format(chrom))
        write_qsub_array_script(script_name, job_name, project, queue, logfile, errfile, array_start, array_end, run_cmd, flag_file)

    #If the file does not already exist try to execute command
    #Define the total outputs
    refs = [os.path.join(out_dir,"ind_{0}_ref.hap.gz".format(i+1)) for i in range(array_start - 1, array_end)]
    already_run = al_done(refs)

    if not already_run:

    #RUN SCRIPT
        try:
            command = "qsub {0}".format(script_name)
            args = shlex.split(command)
            qsub = subprocess.Popen(args, stdout=subprocess.PIPE)
            out = qsub.communicate()[0]
            print (out.strip())

            #CHECK FOR COMPLETION OF TASKS

            #Define flagfiles to look for
            flagfiles = ["{0}.{1}".format(flag_file, i + 1) for i in range(array_start-1, array_end)]
            poll_for_done_mult(flagfiles, 120)
        except KeyboardInterrupt:
            qsub.kill()
        except Exception as e:
            print("[ERROR] Create reference --> Creation of reference files failed. Exiting...]")
            print(e)
            sys.exit()

        #Remove the flagfiles
        for file in flagfiles:
            os.remove(file)
    #Check step variable
    if step == "reference":
        sys.exit()

    ############################################################################
    ############ CREATE VALIDATION AND TARGET ##################################
    ############################################################################

    val_tarlog = os.path.join(log_dir, "val_tarlog")
    if not os.path.exists(val_tarlog):
        os.makedirs(val_tarlog)
    #Define and create output dir for each chromosome and array (validation and test)
    out_val=os.path.join(val_dir, "{0}_chr{1}".format(array, chrom))
    out_target=os.path.join(tar_dir, "{0}_chr{1}".format(array, chrom))

    if not os.path.exists(out_val):
        os.makedirs(out_val)
    if not os.path.exists(out_target):
        os.makedirs(out_target)

    #Define the base name for the outputs file
    val_name=os.path.join(out_val, "ind_${{SGE_TASK_ID}}_{0}_chr{1}".format(array, chrom))
    tar_name=os.path.join(out_target, "ind_${{SGE_TASK_ID}}_{0}_chr{1}".format(array, chrom))

    #Define array positions dir and variables
    array_pos=os.path.join(array_dir, "{0}_{1}_positions.txt".format(array, chrom)) #positions in the array
    array_nopos=os.path.join(array_dir, "{0}_{1}_nopos.txt".format(array, chrom)) #positions not in the array



    #############################
    ## Creating validation set ##
    #############################


    #Subset array positions for individual i
    subset="{0} --gzvcf {2} --indv $id --exclude-positions {1} --recode --recode-INFO-all --stdout | {3} -c > {4}.vcf.gz".format(vcftools, array_pos, vcf_file, bgzip, val_name)

    #Sort
    sort="{0} -c {1}.vcf.gz | {2} -c > {1}_sorted.vcf.gz".format(vcf_sort, val_name, bgzip)

    #Make index
    index="{0} -p vcf {1}_sorted.vcf.gz".format(tabix, val_name)

    #From vcf to impute format
    form="{0} {1} -vcf {2}_sorted.vcf.gz -chr {3} -snps_only -gen {2}_validation".format(perl, val_script, val_name,
                                                                                        chrom)
    #Remove vcfs
    remove="rm {0}.vcf.gz; rm {0}_sorted.vcf.gz*".format(val_name)

    val_cmd = "\n".join([get_id, subset, sort, index, form, remove])

    #########################
    ## Creating target set ##
    #########################

    #Subset array positions for individual i

    subset="{0} --gzvcf {2} --indv $id --positions {1} --recode --recode-INFO-all --stdout  | {3} -c > {4}.vcf.gz".format(vcftools, array_pos, vcf_file, bgzip, tar_name)

    #Sort
    sort="{0} -c {1}.vcf.gz | {2} -c > {1}_sorted.vcf.gz".format(vcf_sort, tar_name, bgzip)

    #Make index
    index="{0} -p vcf {1}_sorted.vcf.gz".format(tabix, tar_name)

    #Tranform to plink to eliminate phasing information
    unphase=unphase_plink(plink, tar_name)

    #Define bed, bim, fam and gentic map files
    bed="{0}_unphased.bed".format(tar_name)
    bim="{0}_unphased.bim".format(tar_name)
    fam="{0}_unphased.fam".format(tar_name)
    in_bed = " ".join([bed,bim,fam])
    refs=" ".join([ref_hap, ref_leg, ref_sample])
    shape_log=os.path.join(val_tarlog, "shapit_log${SGE_TASK_ID}")

    #Phase using shapeit

    phase="{0} --input-bed {1} -M {2} --input-ref {3} -O {4}_target --output-log {5}".format(shapeit, in_bed, gen_map,
                                                                                            refs, tar_name, shape_log)

    #Zip the result

    gzip="gzip {0}_target.haps".format(tar_name)

    #Remove files

    remove="rm {0}_unphased*; rm {0}_sorted*; rm {0}.vcf.gz".format(tar_name)

    tar_cmd = "\n".join([subset, sort, index, unphase, phase, gzip, remove])

    ####################################################################
    ## Joining validation and target commands and writing qsub script ##
    ####################################################################

    run_cmd = "\n".join([val_cmd, tar_cmd])

    #Define the script name and write the script
    script_name=os.path.join(script_dir,"{0}_{1}.sh".format(val_tarname, chrom))
    flag_file = os.path.join(flag_dir, "chr{0}_flagfile".format(chrom))
    #Check if script and out files already exists
    if not os.path.exists(script_name):
        job_name = "{0}_{1}".format(val_tarname, chrom)
        logfile = os.path.join(val_tarlog, "chr{0}_valog.txt".format(chrom))
        errfile = os.path.join(val_tarlog, "chr{0}_valerr.txt".format(chrom))
        write_qsub_array_script(script_name, job_name, project, queue, logfile, errfile, array_start, array_end, run_cmd, flag_file)

    #RUN SCRIPT
    val_files = [os.path.join(out_val,"ind_{0}_{1}_chr{2}_validation.gz".format(i +1, array, chrom)) for i in range(array_start-1, array_end)]
    already_run = al_done(val_files)

    if not already_run:
        try:
            command = "qsub {0}".format(script_name)
            args = shlex.split(command)
            qsub = subprocess.Popen(args, stdout=subprocess.PIPE)
            out = qsub.communicate()[0]
            print (out.strip())

            #CHECK FOR COMPLETION OF TASKS

            #Define flagfiles to look for
            flagfiles = ["{0}.{1}".format(flag_file, i + 1) for i in range(array_start-1,array_end)]
            poll_for_done_mult(flagfiles, 120)
        except KeyboardInterrupt:
            qsub.kill()
        except Exception as e:
            print("[ERROR] Create validation and target --> Creation of val_tar files failed. Exiting...]")
            print(e)
            sys.exit()

        #Remove the flagfiles
        for file in flagfiles:
            os.remove(file)

    ############################################################################
    ############ Write and execute imputation script ###########################
    ############################################################################

    #Create dir for impute logs and outfiles
    imp_log = os.path.join(log_dir, "impute_chr{0}_{1}".format(chrom, array))
    imp_out = os.path.join(imp_res, "chr{0}_{1}".format(chrom, array))

    if not os.path.exists(imp_log):
        os.makedirs(imp_log)
    if not os.path.exists(imp_out):
        os.makedirs(imp_out)

    #Create a directory for every individual

    create_dirs = [os.path.join(imp_out, "indiv{0}".format(i+1)) for i in range(array_start-1, array_end)]
    for path in create_dirs:
        if not os.path.exists(path):
            os.makedirs(path)

    ###################################
    ## Write generic imputation code ##
    ###################################

    start_i="strt=$1"
    right_i="end=$2"
    index="index=$3"
    ind_dir="ind_dir={0}/indiv${{SGE_TASK_ID}}".format(imp_out)
    #Define the file variables

    target_hap = os.path.join(out_target, "ind_${{SGE_TASK_ID}}_{0}_chr{1}_target.haps.gz".format(array, chrom))
    ind_refhap = os.path.join(out_dir, "ind_${SGE_TASK_ID}_ref.hap.gz")
    #previous ind_refleg was created from vcf, new ind_refleg is the same as the unfiltered 1KGP legend.
    #ind_refleg = os.path.join(out_dir, "ind_${SGE_TASK_ID}_ref.legend.gz")
    ind_refleg = ref_leg

    index_file = os.path.join(index_dir, "chr{0}_index.txt".format(chrom))
    out_file = "${ind_dir}/ind_${SGE_TASK_ID}_chunk${index}"

    #Write impute commands
    if imp_mode == "one":
        imp_cmd = impute_oneref(impute2, target_hap, gen_map, ind_refhap, ind_refleg, out_file)
    if imp_mode == "two":
        nat_refhap = os.path.join(second_refdir, "{0}.chr{1}.hap.gz".format(second_basename, chrom))
        nat_refleg = os.path.join(second_refdir, "{0}.chr{1}.legend.gz".format(second_basename, chrom))
        imp_cmd = impute_tworef(impute2, target_hap, gen_map, ind_refhap, ind_refleg, nat_refhap, nat_refleg, out_file)

    run_cmd = "\n".join([start_i, right_i, index, ind_dir, imp_cmd])
    #Write generic impute script
    script_name=os.path.join(script_dir,"{0}_{1}.sh".format(imp_name, chrom))
    flag_file = os.path.join(flag_dir, "flag_chr{0}".format(chrom))
    if not os.path.exists(script_name):
        job_name = "{0}{1}".format(imp_name, chrom)
        logfile = os.path.join(imp_log, "chr{0}_log.txt".format(chrom))
        errfile = os.path.join(imp_log, "chr{0}_log.txt".format(chrom))
        write_qsub_array_script(script_name, job_name, project, queue, logfile, errfile, array_start, array_end, run_cmd, flag_file)

    ##########################################
    ##Submit the multiple imputation scripts##
    ##########################################

    #Read array index
    strt = []
    end = []
    with open(index_file, "r") as file:
        for line in file:
            split=re.split("\t", line)
            if split[0] == str(chrom):
                strt.append(split[1])
                end.append(split[2])

    if len(strt) == len(end):
        index_len = len(strt)
    else:
        print("Sos un boludo")

    #Check if the analysis is already done

    done = []
    k = 0
    for i in range(array_start-1, array_end):
        files = []
        for j in range(index_len):
            files.append(os.path.join(create_dirs[k], "ind_{0}_chunk{1}_summary".format(i+1, j+1)))

        done.append(al_done(files))
        k += 1
    array_len = (array_end - array_start) + 1
    if sum(done) == array_len:
        already_run = True
    else:
        already_run = False

    #Execute imputation scripts

    if not already_run:

        for i in range(0, index_len):
            try:
                command = "qsub {0} {1} {2} {3}".format(script_name, strt[i], end[i], str(i+1))
                args = shlex.split(command)
                qsub = subprocess.Popen(args, stdout=subprocess.PIPE)
                out = qsub.communicate()[0]
                print (out.strip())

                #CHECK FOR COMPLETION OF TASKS

                #Define flagfiles to look for
                flagfiles = ["{0}.{1}".format(flag_file, i + 1) for i in range(array_start -1, array_end)]
                #Poll for done will only check if the particular chunck is done so no need
                #to add another flag
                poll_for_done_mult(flagfiles, 120)
            except KeyboardInterrupt:
                qsub.kill()
            except Exception as e:
                print("[ERROR] Imputation --> Imputation failed. Exiting...]")
                print(e)
                sys.exit()

    #Remove the flagfiles
            for file in flagfiles:
                os.remove(file)

    ###############################################
    ######Scripts para analisis de datos###########
    ###############################################

    #Crete outdir for the summary data, infomaf data and figures.
    #First create the outdir for the chromosome

    sums_dir_chr=os.path.join(sums_dir, "chr{0}_{1}".format(chrom, array))
    infomaf_dir_chr=os.path.join(infomaf_dir, "chr{0}_{1}".format(chrom, array))
    figure_dir_chr=os.path.join(figure_dir, "chr{0}_{1}".format(chrom, array))

    if not os.path.exists(sums_dir_chr):
        os.makedirs(sums_dir_chr)

    #Directories by ancestry
    bind_files = [os.path.join(sums_dir_chr,"ind{0}_completedata.txt".format(i +1)) for i in range(array_start - 1, array_end)]
    already_run = al_done(bind_files)

    #Variables que necesito: imp_out, out_val, array, chrom, ancestry_dir, ind, sums_dir_chr

    shmem = "#$ -pe shmem 10"
    ind = "ind=$SGE_TASK_ID"

    rcom = "{0} {1} --dir={2} --ind=$ind --array={3} --chr={4} --val={5} --out={6}".format(R, bind_script, imp_out, array, chrom, out_val, sums_dir_chr)

    run_cmd = "\n".join([shmem,ind, rcom])

    #Write into a script
    script_name = os.path.join(script_dir, "{0}_{1}.sh".format(bind_name, chrom))
    flag_file = os.path.join(flag_dir, "chr{0}_flagfile".format(chrom))
    if not os.path.exists(script_name):
        job_name = "{0}_{1}".format(bind_name, chrom)
        logfile = os.path.join(log_dir, "chr{0}_bindlog.txt".format(chrom))
        errfile = os.path.join(log_dir, "chr{0}_binderr.txt".format(chrom))
        write_qsub_array_script(script_name, job_name, project, queue, logfile, errfile, array_start, array_end, run_cmd, flag_file)

    #RUN SCRIPT
    #Check if the output dirs are not empty

    if not already_run:
        try:
            command = "qsub {0}".format(script_name)
            args = shlex.split(command)
            qsub = subprocess.Popen(args, stdout=subprocess.PIPE)
            out = qsub.communicate()[0]
            print (out.strip())

            #CHECK FOR COMPLETION OF TASKS

            #Define flagfiles to look for
            flagfiles = ["{0}.{1}".format(flag_file, i + 1) for i in range(array_start-1, array_end)]
            poll_for_done_mult(flagfiles, 120)
        except KeyboardInterrupt:
            qsub.kill()
        except Exception as e:
            print("[ERROR] Create validation and target --> Creation of val_tar files failed. Exiting...]")
            print(e)
            sys.exit()

        #Remove the flagfiles
        for file in flagfiles:
            os.remove(file)

    #After running the binding and

    ##############################################
    ## Script for creating summaries of results ##
    ##############################################

    #The files that I'll need: sums_dir_chr, individual id, sums_dir_chr again

    pops = ["MXL", "CLM", "PUR", "PEL"]
    new_dirs = [os.path.join(sums_dir_chr, pop) for pop in pops]
    for dir in new_dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)

    ind = "ind=$SGE_TASK_ID"

    rcom = "{0} {1} --ind=$ind --dir={2} --out={2}".format(R, sums_script, sums_dir_chr)
    run_cmd = "\n".join([ind, rcom])

    #Write into a script
    script_name = os.path.join(script_dir, "{0}_{1}.sh".format(sums_name, chrom))
    flag_file = os.path.join(flag_dir, "chr{0}_flagfile".format(chrom))
    if not os.path.exists(script_name):
        job_name = "{0}_{1}".format(sums_name, chrom)
        logfile = os.path.join(log_dir, "chr{0}_sumlog.txt".format(chrom))
        errfile = os.path.join(log_dir, "chr{0}_sumerr.txt".format(chrom))
        write_qsub_array_script(script_name, job_name, project, queue, logfile, errfile, array_start, array_end, run_cmd, flag_file)

    #RUN SCRIPT
    #Check if the output dirs are not empty
    already_run = 0

    for dir in new_dirs:
        if len(os.listdir(dir)) != 0:
            already_run = 1


    if not already_run:
        try:
            command = "qsub {0}".format(script_name)
            args = shlex.split(command)
            qsub = subprocess.Popen(args, stdout=subprocess.PIPE)
            out = qsub.communicate()[0]
            print (out.strip())

            #CHECK FOR COMPLETION OF TASKS

            #Define flagfiles to look for
            flagfiles = ["{0}.{1}".format(flag_file, i + 1) for i in range(array_start - 1, array_end)]
            poll_for_done_mult(flagfiles, 120)
        except KeyboardInterrupt:
            qsub.kill()
        except Exception as e:
            print("[ERROR] Create validation and target --> Creation of val_tar files failed. Exiting...]")
            print(e)
            sys.exit()

        #Remove the flagfiles
        for file in flagfiles:
            os.remove(file)

    #################################################
    ## Script for computing concordance of results ##
    #################################################

    #The files that I'll need: sums_dir_chr, individual id, sums_dir_chr again


    ind = "ind=$SGE_TASK_ID"

    rcom = "{0} {1} --ind=$ind --dir={2} --out={2}".format(R, concor_script, sums_dir_chr)
    run_cmd = "\n".join([ind, rcom])

    #Write into a script
    script_name = os.path.join(script_dir, "{0}_{1}.sh".format(concor_name, chrom))
    flag_file = os.path.join(flag_dir, "chr{0}_flagfile".format(chrom))
    if not os.path.exists(script_name):
        job_name = "{0}_{1}".format(concor_name, chrom)
        logfile = os.path.join(log_dir, "chr{0}_concorlog.txt".format(chrom))
        errfile = os.path.join(log_dir, "chr{0}_concorerr.txt".format(chrom))
        write_qsub_array_script(script_name, job_name, project, queue, logfile, errfile, array_start, array_end, run_cmd, flag_file)

    #RUN SCRIPT
    #Check if the output dirs are not empty
    already_run = 0

    #for dir in new_dirs:
    #    if len(os.listdir(dir)) != 0:
    #        already_run = 1


    if not already_run:
        try:
            command = "qsub {0}".format(script_name)
            args = shlex.split(command)
            qsub = subprocess.Popen(args, stdout=subprocess.PIPE)
            out = qsub.communicate()[0]
            print (out.strip())

            #CHECK FOR COMPLETION OF TASKS

            #Define flagfiles to look for
            flagfiles = ["{0}.{1}".format(flag_file, i + 1) for i in range(array_start - 1, array_end)]
            poll_for_done_mult(flagfiles, 120)
        except KeyboardInterrupt:
            qsub.kill()
        except Exception as e:
            print("[ERROR] Computing concordance --> Calculation of concordance failed. Exiting...]")
            print(e)
            sys.exit()

        #Remove the flagfiles
        for file in flagfiles:
            os.remove(file)

    #######################################################
    ## Script for creating info maf summaries of results ##
    #######################################################

    #The files that I'll need: nfomaf_dir_chr, individual id, infomaf_dir_chr again, imp_out, infomaf_name

    pops = ["MXL", "CLM", "PUR", "PEL"]
    new_dirs = [os.path.join(infomaf_dir_chr, pop) for pop in pops]
    for dir in new_dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)

    ind = "ind=$SGE_TASK_ID"

    rcom = "{0} {1} --ind=$ind --dir={2} --out={3}".format(R, infomaf_script, imp_out, infomaf_dir_chr)
    run_cmd = "\n".join([ind, rcom])

    #Write into a script
    script_name = os.path.join(script_dir, "{0}_{1}.sh".format(infomaf_name, chrom))
    flag_file = os.path.join(flag_dir, "chr{0}_flagfile".format(chrom))
    if not os.path.exists(script_name):
        job_name = "{0}_{1}".format(infomaf_name, chrom)
        logfile = os.path.join(log_dir, "chr{0}_infomaflog.txt".format(chrom))
        errfile = os.path.join(log_dir, "chr{0}_infomaferr.txt".format(chrom))
        write_qsub_array_script(script_name, job_name, project, queue, logfile, errfile, array_start, array_end, run_cmd, flag_file)

    #RUN SCRIPT
    #Check if the output dirs are not empty
    already_run = 0

    #Modificar esto mas tarde para que revise bien la condicion
    #for dir in new_dirs:
    #    if len(os.listdir(dir)) != 0:
    #        already_run = 1


    if not already_run:
        try:
            command = "qsub {0}".format(script_name)
            args = shlex.split(command)
            qsub = subprocess.Popen(args, stdout=subprocess.PIPE)
            out = qsub.communicate()[0]
            print (out.strip())

            #CHECK FOR COMPLETION OF TASKS

            #Define flagfiles to look for
            flagfiles = ["{0}.{1}".format(flag_file, i + 1) for i in range(array_start - 1, array_end)]
            poll_for_done_mult(flagfiles, 120)
        except KeyboardInterrupt:
            qsub.kill()
        except Exception as e:
            print("[ERROR] Create validation and target --> Creation of val_tar files failed. Exiting...]")
            print(e)
            sys.exit()

        #Remove the flagfiles
        for file in flagfiles:
            os.remove(file)

    #######################################################
    ## Script for ploting info maf summaries of results ###
    #######################################################

    #Necesito: infomaf_dir_chr, figure_dir_chr

    pops = ["MXL", "CLM", "PUR", "PEL"]
    new_dirs = [os.path.join(figure_dir_chr, pop) for pop in pops]
    for dir in new_dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)


    run_cmd = "{0} {1} --dir={2} --out={3}".format(R, plotmafinfo_script, infomaf_dir_chr, figure_dir_chr)

    #Write into a script
    script_name = os.path.join(script_dir, "{0}_{1}.sh".format(plotmafinfo_name, chrom))
    flag_file = os.path.join(flag_dir, "chr{0}_flagfile".format(chrom))
    if not os.path.exists(script_name):
        job_name = "{0}_{1}".format(plotmafinfo_name, chrom)
        logfile = os.path.join(log_dir, "chr{0}_plotlog.txt".format(chrom))
        errfile = os.path.join(log_dir, "chr{0}_ploterr.txt".format(chrom))
        write_qsub_script(script_name, job_name, project, queue, logfile, errfile, run_cmd, flag_file)

    #RUN SCRIPT
    #Check if the output dirs are not empty
    already_run = 0

    #Modificar esto mas tarde para que revise bien la condicion
    #for dir in new_dirs:
    #    if len(os.listdir(dir)) != 0:
    #        already_run = 1


    if not already_run:
        try:
            command = "qsub {0}".format(script_name)
            args = shlex.split(command)
            qsub = subprocess.Popen(args, stdout=subprocess.PIPE)
            out = qsub.communicate()[0]
            print (out.strip())
            #CHECK FOR COMPLETION OF TASKS
            poll_for_done(flag_file, 120)
        except KeyboardInterrupt:
            qsub.kill()
        except Exception as e:
            print("[ERROR] Create validation and target --> Creation of val_tar files failed. Exiting...]")
            print(e)
            sys.exit()

        #Remove the flagfiles
        os.remove(flag_file)
