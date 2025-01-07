#!/usr/bin/env python3
import os
from os import path
import shutil
import glob
import subprocess
import argparse
import time
from time import sleep
import datetime
import sys
import ntpath
import pandas as pd
import numpy as np
import fileinput
import zlib
import datetime
import sh
import re
from sh import gunzip
import multiprocessing
from pathlib import Path
from Pincho_v01 import *

### master directories
cwd = os.getcwd()
bin_dir = os.path.join(sys.path[0], 'bin')
lib_dir = os.path.join(sys.path[0], 'lib')
pinc_dir = sys.path[0]
conf_file = os.path.join(sys.path[0], 'pincho_v01.conf')

### exporting directories
os.environ["SHANNON_RNA_SEQ_PATH"] = os.path.join(bin_dir, "rcorrector/run_rcorrector.pl")
os.environ["BUSCO_CONFIG_FILE"] = os.path.join(bin_dir, "busco/config/config.ini")
os.environ["TRINITY_HOME"] = os.path.join(bin_dir, "/trinityrnaseq-v2.11.0.FULL/trinityrnaseq-v2.11.0")
os.environ["AUGUSTUS_CONFIG_PATH"] = os.path.join(bin_dir, "Augustus-master/config/")
my_env = os.environ.copy()
os.environ["PATH"] = os.path.join(bin_dir, "RSEM-1.3.1:") + os.path.join(bin_dir, "oases-0.2.09/oases:") + os.path.join(bin_dir, "abyss-2.2.4/bin:") + os.path.join(bin_dir, "Augustus-master/bin:") + os.path.join(bin_dir, "Augustus-master/scripts:") + os.path.join(bin_dir, "velvet_1.2.10:") + (os.path.join(bin_dir, "salmon-1.0.0_linux_x86_64/salmon-latest_linux_x86_64/bin:")) + os.path.join(bin_dir, "blat:") + my_env["PATH"]

### cleaning
trimmomatic_dir = os.path.join(bin_dir, "Trimmomatic-0.39")
rcorrector_dir = os.path.join(bin_dir, "Rcorrector-1.0.4")
cd_hit_dir = os.path.join(bin_dir, "cd-hit-v4.8.1-2019-0228")
transcriptome_tools_dir = os.path.join(bin_dir, "TranscriptomeAssemblyTools-master")

### assemblers
abyss_dir = os.path.join(bin_dir, "abyss-2.2.4")
transabyss_dir = os.path.join(bin_dir, "transabyss-2.0.1")
megahit_dir = os.path.join(bin_dir, "MEGAHIT-1.2.9-Linux-x86_64-static/bin")
oases_script = os.path.join(bin_dir, "oases-0.2.09/oases/scripts/oases_pipeline.py")
shannon_dir = os.path.join(bin_dir, "shannon_cpp-Linux64-v0.4.0")
spades_dir = os.path.join(bin_dir, "SPAdes-3.14.1-Linux")
tadpole_dir = os.path.join(bin_dir, "BBMap_38.86/bbmap")
trinity_dir = os.path.join(bin_dir, "trinityrnaseq-v2.11.0.FULL/trinityrnaseq-v2.11.0")
translig_dir = os.path.join(bin_dir, "TransLiG_1.3")
idba_tran_dir = os.path.join(bin_dir, "idba-1.1.3/bin")
binpacker_dir = os.path.join(bin_dir, "BinPacker_1.0")

### annotation
blast_dir = os.path.join(bin_dir, "ncbi-blast-2.10.0+-x64-linux/ncbi-blast-2.10.0+")
makeblastdb_dir = os.path.join(bin_dir, "ncbi-blast-2.3.0+-x64-linux/ncbi-blast-2.3.0+/bin")

### alingers
bwa_dir = os.path.join(bin_dir, "bwa-0.7.17")
paladin_dir = os.path.join(bin_dir, "paladin-1.4.6")
hisat2_dir = os.path.join(bin_dir, "hisat2-2.1.0-Linux_x86_64/hisat2-2.1.0")
bowtie_bin = 'bowtie2'

### quality assessment
busco_dir = os.path.join(bin_dir, "busco/src/busco")
detonate_dir = os.path.join(bin_dir, "detonate-1.9/rsem-eval")
transrate_dir = os.path.join(bin_dir, "transrate-1.0.3-linux-x86_64")
transrate_read_mapping_dir = os.path.join(bin_dir, "orp-transrate")

### expression analysis
salmon_dir = os.path.join(bin_dir, "salmon-1.0.0_linux_x86_64/salmon-latest_linux_x86_64")
kallisto_dir = os.path.join(bin_dir, "kallisto_linux-v0.46.1/kallisto")

### misc
sratools_dir = os.path.join(bin_dir, "sratoolkit.2.10.8-ubuntu64/bin")
samtools_dir = os.path.join(bin_dir, "samtools-1.10")
mapmitos_dir = os.path.join(bin_dir, "map2mitos")
ancestral_file = os.path.join(lib_dir, "eukaryota_odb10.2019-11-20/eukaryota_odb10/ancestral")
ancestralv_file = os.path.join(lib_dir, "eukaryota_odb10.2019-11-20/eukaryota_odb10/ancestral_variants")
seqkit_dir = os.path.join(bin_dir, "seqkit_linux_amd64")
species_name = "_"

### Load Configuration Page
parser = argparse.ArgumentParser()
parser.add_argument('-conf', metavar='FILE', type=str, required=True, help="pincho_v01 c_fig File")
args = parser.parse_args()

c_fig = {}
with open(args.conf, 'r') as fh_conf:
    for c in fh_conf.readlines():
        c = c.strip()
        if len(c) >0 and c[0] != "#":
            vals = c.split("=")
            c_fig[vals[0].strip()] = vals[1].strip()


### Configurations
evalue          = c_fig['blast_evalue1']
evalue2         = c_fig['blast_evalue2']
evalue3         = c_fig['blast_evalue3']
#sim_runs        = int(c_fig['sim_runs']) #NOTE sim runs must be manually started by user
sim_runs        = 1
memamount       = int(c_fig['nmemory'])
halfmem         = int(int(c_fig['nmemory'])/1)
shannonmem      = int(int(c_fig['nmemory']) * 1000000000)
nfolders        = int(c_fig['nfolders'])
trim_len        = int(c_fig['trim_length'])
megahit_k       = c_fig['megahit_k']
abyss_k         = c_fig['abyss_k']
query_seq       = c_fig['query_seq']
### libraries
busco_lin       = os.path.dirname(c_fig['busco_lin'])
adapter_seqs    = c_fig['adapter_seqs']
annotat_fasta1   = c_fig['annotation_fasta1']
annotat_fasta2   = c_fig['annotation_fasta2']
max_intron      = c_fig['max_intron']
genome_fasta    = c_fig['genome_guide']

def makeblastdb():
    subprocess.run(f"sed -n 2p {check_db} > testing_db.txt", shell=True)
    subprocess.run(f"sed -i 's/[nN]//g;s/[xX]//g;s/[cC]//g;s/[aA]//g;s/[gG]//g;s/[tT]//g;s/ //g;/^$/d' testing_db.txt", shell=True)
    db_out = Path(check_db).stem
    if os.stat("testing_db.txt").st_size == 0:
        subprocess.run(f"{makeblastdb_dir}/makeblastdb -in {check_db} -dbtype nucl -parse_seqids -out {db_out}", shell=True)
    else:
        subprocess.run(f"{makeblastdb_dir}/makeblastdb -in {check_db} -dbtype prot -parse_seqids -out {db_out}", shell=True)
    subprocess.run(f"rm -f testing_db.txt", shell=True)  

if c_fig['blast_db1'] != '0':
    if c_fig['blast_db1'].endswith('.fasta'):
        check_db = c_fig['blast_db1']
        makeblastdb()
        custom_blast_dir     = os.path.dirname(c_fig['blast_db1'])
        blast_name           = Path(c_fig['blast_db1']).stem
    elif c_fig['blast_db1'].endswith('.fa'):
        check_db = c_fig['blast_db1']
        makeblastdb()
        custom_blast_dir     = os.path.dirname(c_fig['blast_db1'])
        blast_name           = Path(c_fig['blast_db1']).stem
    else:
        custom_blast_dir     = os.path.dirname(c_fig['blast_db1'])
        blast_name           = Path(c_fig['blast_db1']).stem

if c_fig['blast_db2'] != '0':
    if c_fig['blast_db2'].endswith('.fasta'):
        check_db = c_fig['blast_db2']
        makeblastdb()
        custom_blast_dir_two     = os.path.dirname(c_fig['blast_db2'])
        blast_name_two           = Path(c_fig['blast_db2']).stem
    elif c_fig['blast_db2'].endswith('.fa'):
        check_db = c_fig['blast_db2']
        makeblastdb()
        custom_blast_dir_two     = os.path.dirname(c_fig['blast_db2'])
        blast_name_two           = Path(c_fig['blast_db2']).stem
    else:
        custom_blast_dir_two     = os.path.dirname(c_fig['blast_db2'])
        blast_name_two           = Path(c_fig['blast_db2']).stem

if c_fig['blast_db3'] != '0':
    if c_fig['blast_db3'].endswith('.fasta'):
        check_db = c_fig['blast_db3']
        makeblastdb()
        custom_blast_dir_three     = os.path.dirname(c_fig['blast_db3'])
        blast_name_three           = Path(c_fig['blast_db3']).stem
    elif c_fig['blast_db3'].endswith('.fa'):
        check_db = c_fig['blast_db3']
        makeblastdb()
        custom_blast_dir_three     = os.path.dirname(c_fig['blast_db3'])
        blast_name_three           = Path(c_fig['blast_db3']).stem
    else:
        custom_blast_dir_three     = os.path.dirname(c_fig['blast_db3'])
        blast_name_three           = Path(c_fig['blast_db3']).stem

if c_fig['abyss_k'] != '0':
    try:
        k0a             = abyss_k.split(',')[0]
    except IndexError:
        k0a = 'none'
    try:
        k1a             = abyss_k.split(',')[1]
    except IndexError:
        k1a = 'none'
    try:
        k2a             = abyss_k.split(',')[2]
    except IndexError:
        k2a = 'none'
    try:
        k3a             = abyss_k.split(',')[3]
    except IndexError:
        k3a = 'none'
    try:
        k4a             = abyss_k.split(',')[4]
    except IndexError:
        k4a = 'none'

spades_k        = c_fig['spades_k']
tadpole_k       = c_fig['tadpole_k']
if c_fig['tadpole_k'] != '0':
    try:
        k0t             = tadpole_k.split(',')[0]
    except IndexError:
        k0t = 'none'
    try:
        k1t             = tadpole_k.split(',')[1]
    except IndexError:
        k1t = 'none'
    try:
        k2t             = tadpole_k.split(',')[2]
    except IndexError:
        k2t = 'none'
    try:        
        k3t             = tadpole_k.split(',')[3]
    except IndexError:
        k3t = 'none'
    try:        
        k4t             = tadpole_k.split(',')[4]
    except IndexError:
        k4t = 'none'
oases_k         = c_fig['oases_k']
if c_fig['oases_k'] != '0':
    k0o             = oases_k.split(',')[0]
    k1o             = oases_k.split(',')[1]
    k2o             = oases_k.split(',')[2]
idba_tran_k     = c_fig['idba_tran_k']
if c_fig['idba_tran_k'] != '0':
    k0i             = idba_tran_k.split(',')[0]
    k1i             = idba_tran_k.split(',')[1]
    k2i             = idba_tran_k.split(',')[2]
rnaspades_k     = c_fig['rnaspades_k']
transabyss_k    = c_fig['transabyss_k']
if c_fig['transabyss_k'] != '0':
    try:
        k0ta            = transabyss_k.split(',')[0]
    except IndexError:
        k0ta = 'none'
    try:
        k1ta            = transabyss_k.split(',')[1]
    except IndexError:
        k1ta = 'none'
    try:
        k2ta            = transabyss_k.split(',')[2]
    except IndexError:
        k2ta = 'none'
    try:
        k3ta            = transabyss_k.split(',')[3]
    except IndexError:
        k3ta = 'none'
    try:
        k4ta            = transabyss_k.split(',')[4]
    except IndexError:
        k4ta = 'none'

if c_fig['blast_type1'] != '0':
    blast_process1 = c_fig['blast_type1']

if c_fig['blast_type2'] != '0':
    blast_process2 = c_fig['blast_type2']

if c_fig['blast_type3'] != '0':
    blast_process3 = c_fig['blast_type3']

### Conditional Configurations
if c_fig['nthreads'] == '0':
    tcoreamount = int(multiprocessing.cpu_count())
else:
    tcoreamount = int(c_fig['nthreads'])
coreamount = int(tcoreamount/sim_runs)
coreamount_abyss = int(coreamount/2)

### Load Commands
def main():
    tfolders = nfolders - 1 #how many folders do you have to run this script on in one directory
    wfolder = -1 #starting number
    autowd = os.getcwd() #save working directory of parent folder

    while(wfolder < tfolders): #while current folder count is less than max number of folders do...
        wfolder = wfolder + 1 #add one to the counter -1 + 1 = 0, where 0 is the first folder in directory
        print(f"{wfolder} / {tfolders}") #display current counter out of max folders
        cawd = os.listdir(autowd)[wfolder] #save name of nth folder to be worked on
        print(cawd) #display name of the nth folder to be worked on
        if path.exists(cawd): #if folder that will be worked on exists...
            os.chdir(cawd) #change directory to the nth folder to be worked on
            print("\n")
            with open(conf_file, "r") as file:
                for line in file:
                    if '=\ty' in line:
                        with open(f'your_pipeline.txt', 'a') as f:
                            f.write((line.split('\t',1)[0])+'()\n')
                            print(line.replace("\n", ""))
                            f.close()
            print("\n")
            manual_reads() if c_fig['read1'] != '0' else sleep(0)
            fastq_dump() if c_fig['SRR'] != '0' else sleep(0)
            trimmomatic() if c_fig['trimmomatic'] == '1' else sleep(0)
            rcorrector() if c_fig['rcorrector'] == '1' else sleep(0)
            megahit() if c_fig['megahit'] == '1' else sleep(0)
            abyss() if c_fig['abyss'] == '1' else sleep(0)
            transabyss() if c_fig['transabyss'] == '1' else sleep(0)
            spades() if c_fig['spades'] == '1' else sleep(0)
            rnaspades() if c_fig['rnaspades'] == '1' else sleep(0)
            tadpole() if c_fig['tadpole'] == '1' else sleep(0)
            shannon() if c_fig['shannon'] == '1' else sleep(0)
            trinity() if c_fig['trinity'] == '1' else sleep(0)
            idba_tran() if c_fig['idba_tran'] == '1' else sleep(0)
            translig() if c_fig['translig'] == '1' else sleep(0)
            binpacker() if c_fig['binpacker'] == '1' else sleep(0)
            oases() if c_fig['oases'] == '1' else sleep(0)            
            transrate() if c_fig['transrate'] == '1' else sleep(0)
            trim_tigs() if c_fig['trim_tigs'] == '1' else sleep(0)
            busco() if c_fig['busco'] == '1' else sleep(0)
            transrate_read_mapping() if c_fig['transrate2'] == '1' else sleep(0)
            hisat2() if c_fig['hisat2'] == '1' else sleep(0)
            cd_hit() if c_fig['cd_hit'] == '1' else sleep(0)
            blast1() if c_fig['blast_type1'] != '0' else sleep(0)
            rsem1() if c_fig['rsem1'] == '1' else sleep(0)
            kallisto1() if c_fig['kallisto1'] == '1' else sleep(0)
            blast2() if c_fig['blast_type2'] != '0' else sleep(0)
            rsem2() if c_fig['rsem2'] == '1' else sleep(0)
            kallisto2() if c_fig['kallisto2'] == '1' else sleep(0)
            blast3() if c_fig['blast_type3'] != '0' else sleep(0)
            cleanup() if c_fig['cleanup'] == '1' else sleep(0)
            os.chdir(autowd)
        else:
            break

if __name__ == '__main__':
    main()

######################################################################################################################

######################################################################################################################
                                            ###   Read Configurations   ###
######################################################################################################################

######################################################################################################################

def r_config():
    global f_reads, r_reads, species_name
    print("r_config...")
    species_name = ntpath.basename(os.path.realpath(os.getcwd()))
    n_files = glob.glob("*.gz")
    latest_file = sorted(n_files, key=os.path.getctime) # sort by creation date
    file1 = latest_file[0] # pull latest before resort
    
    owd = os.getcwd()
    if path.exists(f"{owd}/forward_paired.cor.fq.gz"): #rcor check
        if path.exists(os.getcwd() + "/forward_paired.cor.fq.gz"):
            f_reads = (os.getcwd() + "/forward_paired.cor.fq.gz")
            r_reads = ''
        if path.exists(os.getcwd() + "/reverse_paired.cor.fq.gz"):
            r_reads = (os.getcwd() + "/reverse_paired.cor.fq.gz")
    elif path.exists(f"{owd}/forward_paired.fastq.gz"): #trim check
        if path.exists(os.getcwd() + "/forward_paired.fastq.gz"):
            f_reads = (os.getcwd() + "/forward_paired.fastq.gz")
            r_reads = ''
        if path.exists(os.getcwd() + "/reverse_paired.fastq.gz"):
            r_reads = (os.getcwd() + "/reverse_paired.fastq.gz")
    elif path.exists(f"{owd}/forward.fq.gz"): #trim check
        if path.exists(os.getcwd() + "/forward.fq.gz"):
            f_reads = (os.getcwd() + "/forward.fq.gz")
            r_reads = ''
        if path.exists(os.getcwd() + "/reverse.fq.gz"):
            r_reads = (os.getcwd() + "/reverse.fq.gz")
    else:
        if len(latest_file) == 1:
            file_1_2 = [file1] # list 2 pulled files
            sorted_var = sorted(file_1_2) #sort list alphabetically
            f_reads = sorted_var[0]
            print(f_reads)
            r_reads = ''
        else:
            file2 = latest_file[1] # pull latest before resort
            file_1_2 = [file1,file2] # list 2 pulled files
            sorted_var = sorted(file_1_2) #sort list alphabetically
            f_reads = sorted_var[0]
            print(f_reads)
            r_reads = sorted_var[1]
            print(r_reads)
    adaptive_kmer()

def manual_reads():
    global f_reads, r_reads
    owd = os.getcwd()
    if c_fig['read1'] != '0':
        f_reads = c_fig['read1']
        print(c_fig['read1'])
        shutil.copy2(f_reads, f"{owd}/forward1.gz")
        f_reads = f"{owd}/forward1.gz"
    if c_fig['read2'] != '0':
        r_reads = c_fig['read2']
        print(c_fig['read2'])
        shutil.copy2(r_reads, f"{owd}/reverse2.gz")
        r_reads = f"{owd}/reverse2.gz"

def fastq_dump():
    global f_reads, r_reads
    start_time = time.time()
    print("Performing Fasterq Dump...")
    SRR = c_fig['SRR']
    coreamount_mod = coreamount - 1
    subprocess.run(f"{sratools_dir}/fasterq-dump --split-files -e {coreamount_mod} {SRR}", shell=True)
    subprocess.run(f"sed -i 's/ /_/g' {SRR}*", shell=True)
    subprocess.run(f"pigz {SRR}*", shell=True)
    if path.exists(os.getcwd() + '/' + SRR + ".fastq.gz"):
        f_reads = (os.getcwd() + '/' + SRR + ".fastq.gz")
        r_reads = ''
    if path.exists(os.getcwd() + '/' + SRR + "_1.fastq.gz"):
        f_reads = (os.getcwd() + '/' + SRR + "_1.fastq.gz")
        r_reads = (os.getcwd() + '/' + SRR + "_2.fastq.gz")
    with open('Time_Log.txt', 'a') as f:
        f.write("Fasterq_Dump: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
        f.close()
    r_config()

######################################################################################################################

######################################################################################################################
                                                ###   Kmer Configuration   ###
######################################################################################################################

######################################################################################################################


def adaptive_kmer():
    global q0, q1, q2, q3, q4, qc, qe, qc1, qc2, qc3, qc4, b1
    if not os.path.exists("length.txt"):
        shutil.copy2(f_reads, "forward.gz")
        gunzip("forward.gz")
        subprocess.run("awk '$0 = length; NR==2 { exit }' forward > length.txt", shell=True)
        subprocess.run(f"sed -i 1d length.txt", shell=True)
        subprocess.run(f"rm -f forward", shell=True)    

    a1 = []
    with open('length.txt') as f:
         for line in f:
             data = line.split()
             a1.append(int(data[0]))

    print(*a1)
    b1 = int(*a1)

    if 100 <= int(*a1) <= 1000:
        agg1 = 99
    if 90 <= int(*a1) <= 99:
        agg1 = 89
    if 80 <= int(*a1) <= 89:
        agg1 = 79
    if 70 <= int(*a1) <= 79:
        agg1 = 69
    if 60 <= int(*a1) <= 69:
        agg1 = 59
    if 50 <= int(*a1) <= 59:
        agg1 = 49
    if 40 <= int(*a1) <= 49:
        agg1 = 39
    if 30 <= int(*a1) <= 39:
        agg1 = 29

    q0 = 21
    agg2 = agg1 - q0
    q1 = int(((agg2 + 1)/4) + q0)
    q2 = int(((agg2 + 1)/2) + q0)
    q3 = int(((3*(agg2 + 1))/4) + q0)
    q4 = agg1

    qc = int((q4 - q0)/4)

    if q1 % 2 == 0:
        q1 = q1 - 1
    else:
        q1 = q1
    if q2 % 2 == 0:
        q2 = q2 - 1
    else:
        q2 = q2
    if q3 % 2 == 0:
        q3 = q3 -1
    else:
        q3 = q3

    if qc % 2 == 0:
        qc = qc
    else:
        qc = qc - 1

    qe = int((qc*3) + 21)

    if qe % 2 == 0:
        qe = qe - 1
    else:
        qe = qe

    qc1 = int(q0 + qc)
    qc2 = int(q0 + (2 * qc))
    qc3 = int(q0 + (3 * qc))
    qc4 = int(q0 + (4 * qc))
        

######################################################################################################################

######################################################################################################################
                                                ###   Relics... Not Currently Used   ###
######################################################################################################################

######################################################################################################################
def makeblastdb():
    subprocess.run(f"sed -n 2p {check_db} > testing_db.txt", shell=True)
    subprocess.run(f"sed -i 's/[nN]//g;s/[xX]//g;s/[cC]//g;s/[aA]//g;s/[gG]//g;s/[tT]//g;s/ //g;/^$/d' testing_db.txt", shell=True)
    db_out = Path(check_db).stem
    if os.stat("testing_db.txt").st_size == 0:
        subprocess.run(f"{makeblastdb_dir}/makeblastdb -in {check_db} -dbtype nucl -parse_seqids -out {db_out}", shell=True)
    else:
        subprocess.run(f"{makeblastdb_dir}/makeblastdb -in {check_db} -dbtype prot -parse_seqids -out {db_out}", shell=True)
    subprocess.run(f"rm -f testing_db.txt", shell=True)    

def db_blastn_pull():
    global contig_file, custom_blast_dir, blast_name, species_name, blast_file, db_prot, db_nucl
    species_name = ntpath.basename(os.path.realpath(os.getcwd()))    
    start_time = time.time()
    '''
    add section to do annotation from individual assemblies if merge was not needed (smaller data sets)
    '''
    try: db_nucl
    except NameError: db_nucl = None
    try: db_prot
    except NameError: db_prot = None
    if db_nucl is not None:
        db_out = Path(db_nucl).stem
    if db_prot is not None:
        db_out = Path(db_prot).stem

    if len(glob.glob(f"cd_contig_file")) >= 1:
        contig_file = (os.getcwd() + "/cd_contig_file")
    else:
        contig_file = (os.getcwd() + "/good_merged_assemblies.fasta")

    if not contig_file:
        print("No contig file was found.")
        return #stops program if contig_file is empty/not set

    else:
        print(f"Assembly folder found: running annotation on {contig_file} with {db_out} library")

        wd = os.getcwd()
        os.chdir(f"{lib_dir}/{db_out}")
        filename = ntpath.basename(f"{contig_file}")
        print(filename)

        subprocess.run(f"{blast_dir}/bin/blastn -query {contig_file} -db {db_out} -evalue {evalue} -max_target_seqs 1 -outfmt '7 std qseqid stitle sscinames staxids' -out {species_name}_[{db_out}]{filename}_blast.table -num_threads {coreamount}", shell=True)
        blast_file = (f"{species_name}_[{db_out}]{filename}_blast.table")
        shutil.move(blast_file, wd)
        os.chdir(wd)
    with open('Time_Log.txt', 'a') as f:
        f.write("db_Blastn: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
        f.close()


    owd = os.getcwd()  
    start_time = time.time()

    blast_table_to_fasta()

    with open('Time_Log.txt', 'a') as f:
        f.write("blast_conversion: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
        f.close()

######################################################################################################################

######################################################################################################################
                                                ###   Pre-Processing   ###
######################################################################################################################

######################################################################################################################


# Remove Illumina Adaptors
def trimmomatic(): #input f_reads, r_reads
    global f_reads, r_reads
    r_config()
    if path.exists(os.getcwd() + "/forward_paired.fastq.gz"):
        print("\n\n\nTrimmomatic files already exist\n\n\n ")
    elif not r_reads:
        start_time = time.time()
        shutil.copy2(adapter_seqs, "TruSeq3-SE.fa")
        subprocess.run(f"java -jar {trimmomatic_dir}/trimmomatic-0.39.jar SE -threads {coreamount} -phred33 {f_reads} forward_paired.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36", shell=True)
        print("\n\n\nSE mode\n\n\n ")
        with open('Time_Log.txt', 'a') as f:
            f.write("Trimmomatic: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()    
    else:
        start_time = time.time()
        shutil.copy2(adapter_seqs, "TruSeq3-PE.fa")
        subprocess.run(f"java -jar {trimmomatic_dir}/trimmomatic-0.39.jar PE -threads {coreamount} -phred33 {f_reads} {r_reads} forward_paired.fastq.gz forward_unpaired.fastq.gz reverse_paired.fastq.gz reverse_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36", shell=True)
        print("\n\n\nPE mode\n\n\n ")
        with open('Time_Log.txt', 'a') as f:
            f.write("Trimmomatic: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()    
    if path.exists(os.getcwd() + "/forward_paired.fastq.gz"):
        f_reads = (os.getcwd() + "/forward_paired.fastq.gz")
    if path.exists(os.getcwd() + "/reverse_paired.fastq.gz"):
        r_reads = (os.getcwd() + "/reverse_paired.fastq.gz")
    if path.exists(os.getcwd() + "/forward_paired.fastq"):
        f_reads = (os.getcwd() + "/forward_paired.fastq")
    if path.exists(os.getcwd() + "/reverse_paired.fastq"):
        r_reads = (os.getcwd() + "/reverse_paired.fastq")
    print(f"Forward reads:\n{f_reads}\n")
    print(f"Reverse reads:\n{r_reads}\n")
    subprocess.run(f"rm -f TruSeq3-SE.fa TruSeq3-PE.fa", shell=True)    

# Remove Erroneous Reads
def rcorrector(): #input f_reads, r_reads
    global f_reads, r_reads
    if path.exists(os.getcwd() + "/forward_paired.fastq.gz"):
        sleep(0)
    else:
        r_config()
    if path.exists(os.getcwd() + "/forward_paired.cor.fq.gz"):
        print("\n\n\nRcorrector files already exist\n\n\n ")
    elif not r_reads:
        start_time = time.time()
        subprocess.run(f"{rcorrector_dir}/run_rcorrector.pl -s {f_reads} -t {coreamount}", shell=True)
        with open('Time_Log.txt', 'a') as f:
            f.write("Rcorrector: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
    else:
        start_time = time.time()
        subprocess.run(f"{rcorrector_dir}/run_rcorrector.pl -1 {f_reads} -2 {r_reads} -t {coreamount}", shell=True)
        with open('Time_Log.txt', 'a') as f:
            f.write("Rcorrector: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
    if path.exists(os.getcwd() + "/forward_paired.cor.fq.gz"):
        f_reads = (os.getcwd() + "/forward_paired.cor.fq.gz")
    if path.exists(os.getcwd() + "/reverse_paired.cor.fq.gz"):
        r_reads = (os.getcwd() + "/reverse_paired.cor.fq.gz")
    if path.exists(os.getcwd() + "/forward_paired.cor.fq"):
        f_reads = (os.getcwd() + "/forward_paired.cor.fq")
    if path.exists(os.getcwd() + "/reverse_paired.cor.fq"):
        r_reads = (os.getcwd() + "/reverse_paired.cor.fq")

    print(f"Forward reads:\n{f_reads}\n")
    print(f"Reverse reads:\n{r_reads}\n")


######################################################################################################################

######################################################################################################################
                                                ###   Assemblers   ###
######################################################################################################################

######################################################################################################################
def r_cond(): #Check if Trim or Rcor was run, if not use original reads to drive process
    owd = os.getcwd()
    if path.exists(f"{owd}/forward_paired.fastq.gz"): #trim check
        sleep(0)
    elif path.exists(f"{owd}/forward_paired.cor.fq.gz"): #rcor check
        sleep(0)
    else:
        r_config() #use original reads

def megahit(): #input f_reads, r_reads
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/megahit_contig")) >= 1:
        print("Megahit assembly already present")
    else:
        start_time = time.time()
        if c_fig['megahit_a'] == '1':
            adaptive_kmer()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{megahit_dir}/megahit -t {coreamount} -r {mergedreads} --k-list {q0},{q1},{q2},{q3},{q4} -o megahit_assembly", shell=True)
            elif not r_reads:
                subprocess.run(f"{megahit_dir}/megahit -t {coreamount} -r {f_reads} --k-list {q0},{q1},{q2},{q3},{q4} -o megahit_assembly", shell=True)
            else:
                subprocess.run(f"{megahit_dir}/megahit -t {coreamount} -1 {f_reads} -2 {r_reads} --k-list {q0},{q1},{q2},{q3},{q4} -o megahit_assembly", shell=True)

        elif c_fig['megahit_k'] != '0':
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{megahit_dir}/megahit -t {coreamount} -r {mergedreads} --k-list {megahit_k} -o megahit_assembly", shell=True)
            elif not r_reads:
                subprocess.run(f"{megahit_dir}/megahit -t {coreamount} -r {f_reads} --k-list {megahit_k} -o megahit_assembly", shell=True)
            else:
                subprocess.run(f"{megahit_dir}/megahit -t {coreamount} -1 {f_reads} -2 {r_reads} --k-list {megahit_k} -o megahit_assembly", shell=True)

        else:
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{megahit_dir}/megahit -t {coreamount} -r {mergedreads} -o megahit_assembly", shell=True)
            elif not r_reads:
                subprocess.run(f"{megahit_dir}/megahit -t {coreamount} -r {f_reads} -o megahit_assembly", shell=True)
            else:
                subprocess.run(f"{megahit_dir}/megahit -t {coreamount} -1 {f_reads} -2 {r_reads} -o megahit_assembly", shell=True)

        with open('Time_Log.txt', 'a') as f:
            f.write(f"Megahit: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/megahit_assembly/final.contigs.fa", "Assemblies_Folder/megahit_contig")
        shutil.rmtree(f"{owd}/megahit_assembly")


def abyss(): #input f_reads, r_reads
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/abyss_contig")) >= 1:
        print("ABySS assembly 1 already present.")
    else:
        start_time = time.time()
        if c_fig['abyss_a'] == '1':
            adaptive_kmer()
            ########################################RUN1##############################################################################
            os.mkdir(f"{owd}/abyss_assembly1")
            os.chdir(f"{owd}/abyss_assembly1")
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q0} name=abyss_run se='{mergedreads}'", shell=True)
            elif not r_reads:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q0} name=abyss_run se='{f_reads}'", shell=True)
            else:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q0} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
            os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
            os.chdir(owd)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"ABySS_{q0}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            if not os.path.exists(f"{owd}/Assemblies_Folder"):
                os.makedirs(f"{owd}/Assemblies_Folder")
            shutil.copy2(f"{owd}/abyss_assembly1/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs1")
            shutil.rmtree(f"{owd}/abyss_assembly1")
            ########################################RUN2##############################################################################
            os.mkdir(f"{owd}/abyss_assembly2")
            os.chdir(f"{owd}/abyss_assembly2")
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q1} name=abyss_run se='{mergedreads}'", shell=True)
            elif not r_reads:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q1} name=abyss_run se='{f_reads}'", shell=True)
            else:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q1} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
            os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
            os.chdir(owd)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"ABySS_{q1}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            shutil.copy2(f"{owd}/abyss_assembly2/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs2")
            shutil.rmtree(f"{owd}/abyss_assembly2")
            ########################################RUN3##############################################################################
            os.mkdir(f"{owd}/abyss_assembly3")
            os.chdir(f"{owd}/abyss_assembly3")
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q2} name=abyss_run se='{mergedreads}'", shell=True)
            elif not r_reads:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q2} name=abyss_run se='{f_reads}'", shell=True)
            else:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q2} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
            os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
            os.chdir(owd)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"ABySS_{q2}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            shutil.copy2(f"{owd}/abyss_assembly3/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs3")
            shutil.rmtree(f"{owd}/abyss_assembly3")
            ########################################RUN4##############################################################################
            os.mkdir(f"{owd}/abyss_assembly4")
            os.chdir(f"{owd}/abyss_assembly4")
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q3} name=abyss_run se='{mergedreads}'", shell=True)
            elif not r_reads:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q3} name=abyss_run se='{f_reads}'", shell=True)
            else:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q3} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
            os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
            os.chdir(owd)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"ABySS_{q3}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            shutil.copy2(f"{owd}/abyss_assembly4/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs4")
            shutil.rmtree(f"{owd}/abyss_assembly4")
            ########################################RUN5##############################################################################
            os.mkdir(f"{owd}/abyss_assembly5")
            os.chdir(f"{owd}/abyss_assembly5")
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q4} name=abyss_run se='{mergedreads}'", shell=True)
            elif not r_reads:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q4} name=abyss_run se='{f_reads}'", shell=True)
            else:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={q4} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
            os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
            os.chdir(owd)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"ABySS_{q4} %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            shutil.copy2(f"{owd}/abyss_assembly5/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs5")
            shutil.rmtree(f"{owd}/abyss_assembly5")

            os.chdir(f"{owd}/Assemblies_Folder")
            subprocess.run(f"cat abyss_contigs* > abyss_contigs6.fasta", shell=True)
            subprocess.run(f"{tadpole_dir}/sortbyname.sh in=abyss_contigs6.fasta out=abyss_contigs7.fasta length descending", shell=True)
            shutil.copy2(f"{owd}/Assemblies_Folder/abyss_contigs7.fasta", f"{owd}/Assemblies_Folder/abyss_contig")
            subprocess.run(f"rm -f abyss_contigs*", shell=True)
            os.chdir(owd)

        elif c_fig['abyss_k'] != '0':
            ########################################RUN1##############################################################################
            os.mkdir(f"{owd}/abyss_assembly1")
            os.chdir(f"{owd}/abyss_assembly1")
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k0a} name=abyss_run se='{mergedreads}'", shell=True)
            elif not r_reads:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k0a} name=abyss_run se='{f_reads}'", shell=True)
            else:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k0a} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
            os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
            os.chdir(owd)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"ABySS_{k0a}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            if not os.path.exists(f"{owd}/Assemblies_Folder"):
                os.makedirs(f"{owd}/Assemblies_Folder")
            shutil.copy2(f"{owd}/abyss_assembly1/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs1")
            shutil.rmtree(f"{owd}/abyss_assembly1")
            ########################################RUN2##############################################################################
            if k1a != 'none':
                os.mkdir(f"{owd}/abyss_assembly2")
                os.chdir(f"{owd}/abyss_assembly2")
                start_time = time.time()
                try: mergedreads
                except NameError: mergedreads = None
                if mergedreads is not None:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k1a} name=abyss_run se='{mergedreads}'", shell=True)
                elif not r_reads:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k1a} name=abyss_run se='{f_reads}'", shell=True)
                else:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k1a} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
                os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
                os.chdir(owd)
                with open('Time_Log.txt', 'a') as f:
                    f.write(f"ABySS_{k1a}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                    f.close()
                shutil.copy2(f"{owd}/abyss_assembly2/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs2")
                shutil.rmtree(f"{owd}/abyss_assembly2")
            ########################################RUN3##############################################################################
            if k2a != 'none':
                os.mkdir(f"{owd}/abyss_assembly3")
                os.chdir(f"{owd}/abyss_assembly3")
                start_time = time.time()
                try: mergedreads
                except NameError: mergedreads = None
                if mergedreads is not None:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k2a} name=abyss_run se='{mergedreads}'", shell=True)
                elif not r_reads:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k2a} name=abyss_run se='{f_reads}'", shell=True)
                else:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k2a} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
                os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
                os.chdir(owd)
                with open('Time_Log.txt', 'a') as f:
                    f.write(f"ABySS_{k2a}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                    f.close()
                shutil.copy2(f"{owd}/abyss_assembly3/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs3")
                shutil.rmtree(f"{owd}/abyss_assembly3")
            ########################################RUN4##############################################################################
            if k3a != 'none':
                os.mkdir(f"{owd}/abyss_assembly4")
                os.chdir(f"{owd}/abyss_assembly4")
                start_time = time.time()
                try: mergedreads
                except NameError: mergedreads = None
                if mergedreads is not None:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k3a} name=abyss_run se='{mergedreads}'", shell=True)
                elif not r_reads:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k3a} name=abyss_run se='{f_reads}'", shell=True)
                else:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k3a} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
                os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
                os.chdir(owd)
                with open('Time_Log.txt', 'a') as f:
                    f.write(f"ABySS_{k3a}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                    f.close()
                shutil.copy2(f"{owd}/abyss_assembly4/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs4")
                shutil.rmtree(f"{owd}/abyss_assembly4")
            ########################################RUN5##############################################################################
            if k4a != 'none':
                os.mkdir(f"{owd}/abyss_assembly5")
                os.chdir(f"{owd}/abyss_assembly5")
                start_time = time.time()
                try: mergedreads
                except NameError: mergedreads = None
                if mergedreads is not None:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k4a} name=abyss_run se='{mergedreads}'", shell=True)
                elif not r_reads:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k4a} name=abyss_run se='{f_reads}'", shell=True)
                else:
                    subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k={k4a} name=abyss_run in='{f_reads} {r_reads}'", shell=True)
                os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
                os.chdir(owd)
                with open('Time_Log.txt', 'a') as f:
                    f.write(f"ABySS_{k4a} %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                    f.close()
                shutil.copy2(f"{owd}/abyss_assembly5/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contigs5")
                shutil.rmtree(f"{owd}/abyss_assembly5")

            os.chdir(f"{owd}/Assemblies_Folder")
            subprocess.run(f"cat abyss_contigs* > abyss_contigs6.fasta", shell=True)
            subprocess.run(f"{tadpole_dir}/sortbyname.sh in=abyss_contigs6.fasta out=abyss_contigs7.fasta length descending", shell=True)
            shutil.copy2(f"{owd}/Assemblies_Folder/abyss_contigs7.fasta", f"{owd}/Assemblies_Folder/abyss_contig")
            subprocess.run(f"rm -f abyss_contigs*", shell=True)
            os.chdir(owd)

        else:
            ########################################RUN1##############################################################################
            os.mkdir(f"{owd}/abyss_assembly1")
            os.chdir(f"{owd}/abyss_assembly1")
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k=23 name=abyss_run se='{mergedreads}'", shell=True)
            elif not r_reads:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k=23 name=abyss_run se='{f_reads}'", shell=True)
            else:
                subprocess.run(f"{abyss_dir}/bin/abyss-pe j={coreamount} k=23 name=abyss_run in='{f_reads} {r_reads}'", shell=True)
            os.rename('abyss_run-unitigs.fa','abyss_run-contigs.fa')
            os.chdir(owd)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"ABySS: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            if not os.path.exists(f"{owd}/Assemblies_Folder"):
                os.makedirs(f"{owd}/Assemblies_Folder")
            shutil.copy2(f"{owd}/abyss_assembly1/abyss_run-contigs.fa", "Assemblies_Folder/abyss_contig")
            shutil.rmtree(f"{owd}/abyss_assembly1")


def transabyss(): #input f_reads, r_reads
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/transabyss_contig")) >= 1:
        print("Trans-ABySS assembly 1 already present")
    else:
        if c_fig['transabyss_a'] == '1':
            adaptive_kmer()
            ########################################RUN1##############################################################################
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q0} --outdir transabyss_assembly1 --name transabyss_contig --se {mergedreads}", shell=True)
            elif not r_reads:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q0} --outdir transabyss_assembly1 --name transabyss_contig --se {f_reads}", shell=True)
            else:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q0} --outdir transabyss_assembly1 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"Trans-ABySS_{q0}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            if not os.path.exists(f"{owd}/Assemblies_Folder"):
                os.makedirs(f"{owd}/Assemblies_Folder")
            shutil.copy2(f"{owd}/transabyss_assembly1/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs1")
            shutil.rmtree(f"{owd}/transabyss_assembly1")
            ########################################RUN2##############################################################################
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q1} --outdir transabyss_assembly2 --name transabyss_contig --se {mergedreads}", shell=True)
            elif not r_reads:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q1} --outdir transabyss_assembly2 --name transabyss_contig --se {f_reads}", shell=True)
            else:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q1} --outdir transabyss_assembly2 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"Trans-ABySS_{q1}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            shutil.copy2(f"{owd}/transabyss_assembly2/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs2")
            shutil.rmtree(f"{owd}/transabyss_assembly2")
            ########################################RUN3##############################################################################
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q2} --outdir transabyss_assembly3 --name transabyss_contig --se {mergedreads}", shell=True)
            elif not r_reads:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q2} --outdir transabyss_assembly3 --name transabyss_contig --se {f_reads}", shell=True)
            else:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q2} --outdir transabyss_assembly3 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"Trans-ABySS_{q2}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            shutil.copy2(f"{owd}/transabyss_assembly3/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs3")
            shutil.rmtree(f"{owd}/transabyss_assembly3")
            ########################################RUN4##############################################################################
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q3} --outdir transabyss_assembly4 --name transabyss_contig --se {mergedreads}", shell=True)
            elif not r_reads:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q3} --outdir transabyss_assembly4 --name transabyss_contig --se {f_reads}", shell=True)
            else:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q3} --outdir transabyss_assembly4 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"Trans-ABySS_{q3}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            shutil.copy2(f"{owd}/transabyss_assembly4/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs4")
            shutil.rmtree(f"{owd}/transabyss_assembly4")
            ########################################RUN5##############################################################################
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q4} --outdir transabyss_assembly5 --name transabyss_contig --se {mergedreads}", shell=True)
            elif not r_reads:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q4} --outdir transabyss_assembly5 --name transabyss_contig --se {f_reads}", shell=True)
            else:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {q4} --outdir transabyss_assembly5 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"Trans-ABySS_{q4}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            shutil.copy2(f"{owd}/transabyss_assembly5/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs5")
            shutil.rmtree(f"{owd}/transabyss_assembly5")
            
            os.chdir(f"{owd}/Assemblies_Folder")
            subprocess.run(f"{transabyss_dir}/transabyss-merge --threads {coreamount} --mink {q0} --maxk {q4} --prefixes k{q0} k{q1} k{q2} k{q3} k{q4} --out transabyss_contigs6.fasta transabyss_contigs1 transabyss_contigs2 transabyss_contigs3 transabyss_contigs4 transabyss_contigs5", shell=True)        
            #subprocess.run(f"cat transabyss_contigs* > transabyss_contigs6.fasta", shell=True)
            subprocess.run(f"{tadpole_dir}/sortbyname.sh in=transabyss_contigs6.fasta out=transabyss_contigs7.fasta length descending", shell=True)
            shutil.copy2(f"{owd}/Assemblies_Folder/transabyss_contigs7.fasta", "transabyss_contig")
            subprocess.run(f"rm -f transabyss_contigs*", shell=True)
            os.chdir(owd)
            subprocess.run(f"sed -i 's/ /_length_/' {owd}/Assemblies_Folder/transabyss_contig", shell=True)


        elif c_fig['transabyss_k'] != '0':
            ########################################RUN1##############################################################################
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k0ta} --outdir transabyss_assembly1 --name transabyss_contig --se {mergedreads}", shell=True)
            elif not r_reads:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k0ta} --outdir transabyss_assembly1 --name transabyss_contig --se {f_reads}", shell=True)
            else:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k0ta} --outdir transabyss_assembly1 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"Trans-ABySS_{k0ta}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            if not os.path.exists(f"{owd}/Assemblies_Folder"):
                os.makedirs(f"{owd}/Assemblies_Folder")
            shutil.copy2(f"{owd}/transabyss_assembly1/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs1")
            shutil.rmtree(f"{owd}/transabyss_assembly1")
            ########################################RUN2##############################################################################
            if k1ta != 'none':
                start_time = time.time()
                try: mergedreads
                except NameError: mergedreads = None
                if mergedreads is not None:
                    subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k1ta} --outdir transabyss_assembly2 --name transabyss_contig --se {mergedreads}", shell=True)
                elif not r_reads:
                    subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k1ta} --outdir transabyss_assembly2 --name transabyss_contig --se {f_reads}", shell=True)
                else:
                    subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k1ta} --outdir transabyss_assembly2 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
                with open('Time_Log.txt', 'a') as f:
                    f.write(f"Trans-ABySS_{k1ta}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                    f.close()
                shutil.copy2(f"{owd}/transabyss_assembly2/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs2")
                shutil.rmtree(f"{owd}/transabyss_assembly2")
                if k2ta != 'none':
                    start_time = time.time()
                    try: mergedreads
                    except NameError: mergedreads = None
                    if mergedreads is not None:
                        subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k2ta} --outdir transabyss_assembly3 --name transabyss_contig --se {mergedreads}", shell=True)
                    elif not r_reads:
                        subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k2ta} --outdir transabyss_assembly3 --name transabyss_contig --se {f_reads}", shell=True)
                    else:
                        subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k2ta} --outdir transabyss_assembly3 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
                    with open('Time_Log.txt', 'a') as f:
                        f.write(f"Trans-ABySS_{k2ta}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                        f.close()
                    shutil.copy2(f"{owd}/transabyss_assembly3/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs3")
                    shutil.rmtree(f"{owd}/transabyss_assembly3")
                    if k3ta != 'none':
                        start_time = time.time()
                        try: mergedreads
                        except NameError: mergedreads = None
                        if mergedreads is not None:
                            subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k3ta} --outdir transabyss_assembly4 --name transabyss_contig --se {mergedreads}", shell=True)
                        elif not r_reads:
                            subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k3ta} --outdir transabyss_assembly4 --name transabyss_contig --se {f_reads}", shell=True)
                        else:
                            subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k3ta} --outdir transabyss_assembly4 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
                        with open('Time_Log.txt', 'a') as f:
                            f.write(f"Trans-ABySS_{k3ta}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                            f.close()
                        shutil.copy2(f"{owd}/transabyss_assembly4/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs4")
                        shutil.rmtree(f"{owd}/transabyss_assembly4")
                        if k4ta != 'none':
                            start_time = time.time()
                            try: mergedreads
                            except NameError: mergedreads = None
                            if mergedreads is not None:
                                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k4ta} --outdir transabyss_assembly5 --name transabyss_contig --se {mergedreads}", shell=True)
                            elif not r_reads:
                                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k4ta} --outdir transabyss_assembly5 --name transabyss_contig --se {f_reads}", shell=True)
                            else:
                                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} -k {k4ta} --outdir transabyss_assembly5 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
                            with open('Time_Log.txt', 'a') as f:
                                f.write(f"Trans-ABySS_{k4ta}: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                                f.close()
                            shutil.copy2(f"{owd}/transabyss_assembly5/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contigs5")
                            shutil.rmtree(f"{owd}/transabyss_assembly5")

                            os.chdir(f"{owd}/Assemblies_Folder")
                            subprocess.run(f"{transabyss_dir}/transabyss-merge --threads {coreamount} --mink {k0ta} --maxk {k4ta} --prefixes k{k0ta} k{k1ta} k{k2ta} k{k3ta} k{k4ta} --out transabyss_contigs6.fasta transabyss_contigs1 transabyss_contigs2 transabyss_contigs3 transabyss_contigs4 transabyss_contigs5", shell=True)        
                            subprocess.run(f"{tadpole_dir}/sortbyname.sh in=transabyss_contigs6.fasta out=transabyss_contigs7.fasta length descending", shell=True)
                            shutil.copy2(f"{owd}/Assemblies_Folder/transabyss_contigs7.fasta", "transabyss_contig")
                            subprocess.run(f"rm -f transabyss_contigs*", shell=True)
                            os.chdir(owd)
                        else:
                            os.chdir(f"{owd}/Assemblies_Folder")
                            subprocess.run(f"{transabyss_dir}/transabyss-merge --threads {coreamount} --mink {k0ta} --maxk {k3ta} --prefixes k{k0ta} k{k1ta} k{k2ta} k{k3ta} --out transabyss_contigs6.fasta transabyss_contigs1 transabyss_contigs2 transabyss_contigs3 transabyss_contigs4", shell=True)        
                            subprocess.run(f"{tadpole_dir}/sortbyname.sh in=transabyss_contigs6.fasta out=transabyss_contigs7.fasta length descending", shell=True)
                            shutil.copy2(f"{owd}/Assemblies_Folder/transabyss_contigs7.fasta", "transabyss_contig")
                            subprocess.run(f"rm -f transabyss_contigs*", shell=True)
                            os.chdir(owd)
                    else:
                        os.chdir(f"{owd}/Assemblies_Folder")
                        subprocess.run(f"{transabyss_dir}/transabyss-merge --threads {coreamount} --mink {k0ta} --maxk {k2ta} --prefixes k{k0ta} k{k1ta} k{k2ta} --out transabyss_contigs6.fasta transabyss_contigs1 transabyss_contigs2 transabyss_contigs3", shell=True)        
                        subprocess.run(f"{tadpole_dir}/sortbyname.sh in=transabyss_contigs6.fasta out=transabyss_contigs7.fasta length descending", shell=True)
                        shutil.copy2(f"{owd}/Assemblies_Folder/transabyss_contigs7.fasta", "transabyss_contig")
                        subprocess.run(f"rm -f transabyss_contigs*", shell=True)
                        os.chdir(owd)
                else:
                    os.chdir(f"{owd}/Assemblies_Folder")
                    subprocess.run(f"{transabyss_dir}/transabyss-merge --threads {coreamount} --mink {k0ta} --maxk {k1ta} --prefixes k{k0ta} k{k1ta} --out transabyss_contigs6.fasta transabyss_contigs1 transabyss_contigs2", shell=True)        
                    subprocess.run(f"{tadpole_dir}/sortbyname.sh in=transabyss_contigs6.fasta out=transabyss_contigs7.fasta length descending", shell=True)
                    shutil.copy2(f"{owd}/Assemblies_Folder/transabyss_contigs7.fasta", "transabyss_contig")
                    subprocess.run(f"rm -f transabyss_contigs*", shell=True)
                    os.chdir(owd)
            else:
                os.chdir(f"{owd}/Assemblies_Folder")
                subprocess.run(f"cat transabyss_contigs* > transabyss_contigs6.fasta", shell=True)
                subprocess.run(f"{tadpole_dir}/sortbyname.sh in=transabyss_contigs6.fasta out=transabyss_contigs7.fasta length descending", shell=True)
                shutil.copy2(f"{owd}/Assemblies_Folder/transabyss_contigs7.fasta", "transabyss_contig")
                subprocess.run(f"rm -f transabyss_contigs*", shell=True)
                os.chdir(owd)
            ########################################RUN3##############################################################################

            ########################################RUN4##############################################################################

            ########################################RUN5##############################################################################

            subprocess.run(f"sed -i 's/ /_length_/' {owd}/Assemblies_Folder/transabyss_contig", shell=True)

        else:
            start_time = time.time()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} --outdir transabyss_assembly1 --name transabyss_contig --se {mergedreads}", shell=True)
            elif not r_reads:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} --outdir transabyss_assembly1 --name transabyss_contig --se {f_reads}", shell=True)
            else:
                subprocess.run(f"{transabyss_dir}/transabyss --threads {coreamount} --outdir transabyss_assembly1 --name transabyss_contig --pe {f_reads} {r_reads}", shell=True)
            with open('Time_Log.txt', 'a') as f:
                f.write(f"Trans-ABySS: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
                f.close()
            if not os.path.exists(f"{owd}/Assemblies_Folder"):
                os.makedirs(f"{owd}/Assemblies_Folder")
            shutil.copy2(f"{owd}/transabyss_assembly1/transabyss_contig-final.fa", "Assemblies_Folder/transabyss_contig")
            shutil.rmtree(f"{owd}/transabyss_assembly1")
            subprocess.run(f"sed -i 's/ /_length_/' {owd}/Assemblies_Folder/transabyss_contig", shell=True)



def spades(): #input f_reads, r_reads
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/spades_contig")) >= 1:
        print("SPAdes assembly already present")
    else:
        start_time = time.time()
        if c_fig['spades_a'] == '1':
            adaptive_kmer()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{spades_dir}/bin/spades.py --only-assembler -k {q0},{q1},{q2},{q3},{q4} -s {mergedreads} -t {coreamount} -o spades_assembly", shell=True)
            elif not r_reads:
                subprocess.run(f"{spades_dir}/bin/spades.py --only-assembler -k {q0},{q1},{q2},{q3},{q4} -s {f_reads} -t {coreamount} -o spades_assembly", shell=True)
            else:
                subprocess.run(f"{spades_dir}/bin/spades.py --only-assembler -k {q0},{q1},{q2},{q3},{q4} -1 {f_reads} -2 {r_reads} -t {coreamount} -o spades_assembly", shell=True)

        elif c_fig['spades_k'] != '0':
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{spades_dir}/bin/spades.py --only-assembler -k {spades_k} -s {mergedreads} -t {coreamount} -o spades_assembly", shell=True)
            elif not r_reads:
                subprocess.run(f"{spades_dir}/bin/spades.py --only-assembler -k {spades_k} -s {f_reads} -t {coreamount} -o spades_assembly", shell=True)
            else:
                subprocess.run(f"{spades_dir}/bin/spades.py --only-assembler -k {spades_k} -1 {f_reads} -2 {r_reads} -t {coreamount} -o spades_assembly", shell=True)

        else:
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{spades_dir}/bin/spades.py --only-assembler -s {mergedreads} -t {coreamount} -o spades_assembly", shell=True)
            elif not r_reads:
                subprocess.run(f"{spades_dir}/bin/spades.py --only-assembler -s {f_reads} -t {coreamount} -o spades_assembly", shell=True)
            else:
                subprocess.run(f"{spades_dir}/bin/spades.py --only-assembler -1 {f_reads} -2 {r_reads} -t {coreamount} -o spades_assembly", shell=True)

        with open('Time_Log.txt', 'a') as f:
            f.write(f"SPAdes: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/spades_assembly/contigs.fasta", "Assemblies_Folder/spades_contig")
        shutil.rmtree(f"{owd}/spades_assembly")


def rnaspades(): #input f_reads, r_reads
    owd = os.getcwd()
    r_cond()
    q4_spades = q4-2
    if len(glob.glob("Assemblies_Folder/rnaspades_contig")) >= 1:
        print("rnaSPAdes assembly already present")
    else:
        start_time = time.time()       
        if c_fig['rnaspades_a'] == '1':
            adaptive_kmer()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{spades_dir}/bin/rnaspades.py --only-assembler -k {q0},{q1},{q2},{q3},{q4_spades} -s {mergedreads} -t {coreamount} -m {memamount} -o rnaspades_assembly", shell=True)
            elif not r_reads:
                subprocess.run(f"{spades_dir}/bin/rnaspades.py --only-assembler -k {q0},{q1},{q2},{q3},{q4_spades} -s {f_reads} -t {coreamount} -m {memamount} -o rnaspades_assembly", shell=True)
            else:
                subprocess.run(f"{spades_dir}/bin/rnaspades.py --only-assembler -k {q0},{q1},{q2},{q3},{q4_spades} -1 {f_reads} -2 {r_reads} -t {coreamount} -m {memamount} -o rnaspades_assembly", shell=True)

        elif c_fig['rnaspades_k'] != '0':
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{spades_dir}/bin/rnaspades.py --only-assembler -k {rnaspades_k} -s {mergedreads} -t {coreamount} -m {memamount} -o rnaspades_assembly", shell=True)
            elif not r_reads:
                subprocess.run(f"{spades_dir}/bin/rnaspades.py --only-assembler -k {rnaspades_k} -s {f_reads} -t {coreamount} -m {memamount} -o rnaspades_assembly", shell=True)
            else:
                subprocess.run(f"{spades_dir}/bin/rnaspades.py --only-assembler -k {rnaspades_k} -1 {f_reads} -2 {r_reads} -t {coreamount} -m {memamount} -o rnaspades_assembly", shell=True)

        else:
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{spades_dir}/bin/rnaspades.py --only-assembler -k 23 -s {mergedreads} -t {coreamount} -m {memamount} -o rnaspades_assembly", shell=True)
            elif not r_reads:
                subprocess.run(f"{spades_dir}/bin/rnaspades.py --only-assembler -k 23 -s {f_reads} -t {coreamount} -m {memamount} -o rnaspades_assembly", shell=True)
            else:
                subprocess.run(f"{spades_dir}/bin/rnaspades.py --only-assembler -k 23 -1 {f_reads} -2 {r_reads} -t {coreamount} -m {memamount} -o rnaspades_assembly", shell=True)


        with open('Time_Log.txt', 'a') as f:
            f.write(f"rnaSPAdes: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/rnaspades_assembly/transcripts.fasta", f"{owd}/Assemblies_Folder/rnaspades_contig")
        shutil.rmtree(f"{owd}/rnaspades_assembly")


def oases(): #input f_reads, r_reads
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/oases_contig")) >= 1:
        print("Oases assembly already present")
    else:
        print("Oases in Progress...")
        start_time = time.time()
        if c_fig['oases_a'] == '1':
            adaptive_kmer()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{oases_script} -m {q0} -M {qe} -s {qc} -o oases_assembly -d ' -fastq.gz {mergedreads} -interleaved'", shell=True)
            elif not r_reads:
                subprocess.run(f"{oases_script} -m {q0} -M {qe} -s {qc} -o oases_assembly -d ' -fastq.gz {f_reads} '", shell=True)
            else:
                if b1 < 124:
                    subprocess.run(f"{oases_script} -m {q0} -M {qe} -s {qc} -o oases_assembly -d ' -short -fastq.gz {f_reads} {r_reads} -separate '", shell=True)
                else:
                    subprocess.run(f"{oases_script} -m {q0} -M {qe} -s {qc} -o oases_assembly -d ' -long -fastq.gz {f_reads} {r_reads} -separate '", shell=True)

        elif c_fig['oases_k'] != '0':
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{oases_script} -m {k0o} -M {k2o} -s {k1o} -o oases_assembly -d ' -fastq.gz {mergedreads} -interleaved'", shell=True)
            elif not r_reads:
                subprocess.run(f"{oases_script} -m {k0o} -M {k2o} -s {k1o} -o oases_assembly -d ' -fastq.gz {f_reads} '", shell=True)
            else:
                if b1 < 124:
                    subprocess.run(f"{oases_script} -m {k0o} -M {k2o} -s {k1o} -o oases_assembly -d ' -short -fastq.gz {f_reads} {r_reads} -separate '", shell=True)
                else:
                    subprocess.run(f"{oases_script} -m {k0o} -M {k2o} -s {k1o} -o oases_assembly -d ' -long -fastq.gz {f_reads} {r_reads} -separate '", shell=True)

        else:
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{oases_script} -o oases_assembly -d ' -fastq.gz {mergedreads} -interleaved'", shell=True)
            elif not r_reads:
                subprocess.run(f"{oases_script} -o oases_assembly -d ' -fastq.gz {f_reads} '", shell=True)
            else:
                if b1 < 124:
                    subprocess.run(f"{oases_script} -o oases_assembly -d ' -short -fastq.gz {f_reads} {r_reads} -separate '", shell=True)
                else:
                    subprocess.run(f"{oases_script} -o oases_assembly -d ' -long -fastq.gz {f_reads} {r_reads} -separate '", shell=True)

        with open('Time_Log.txt', 'a') as f:
            f.write(f"Oases: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        for oases_assembly in glob.glob(f"{owd}/oases_assembly_*"):
            shutil.move(oases_assembly, f"{owd}/oases_assemblyMerged") #move all oases directories
        print("Completed Oases Script")
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/oases_assemblyMerged/transcripts.fa", "Assemblies_Folder/oases_contig")
        shutil.rmtree(f"{owd}/oases_assemblyMerged")


def trinity(): #input f_reads, r_reads
    owd = os.getcwd()
    r_cond()
    print(genome_fasta)
    if len(glob.glob("Assemblies_Folder/trinity_contig")) >= 1:
        print("Trinity assembly already present")
    elif genome_fasta != "0" and "genome fasta":
        start_time = time.time()
        shutil.copy2(genome_fasta, "hisat_alignment_ref")
        subprocess.run(f"{hisat2_dir}/hisat2-build hisat_alignment_ref hisat_index", shell=True)
        if not r_reads:
            subprocess.run(f"{hisat2_dir}/hisat2 -x hisat_index -U {f_reads} -p {coreamount} -S aligned_reads.sam --summary-file hisat2_summary --al-conc aligned_reads", shell=True)
        else:
            subprocess.run(f"{hisat2_dir}/hisat2 -x hisat_index -1 {f_reads} -2 {r_reads} -p {coreamount} -S aligned_reads.sam --summary-file hisat2_summary --al-conc aligned_reads", shell=True)
        subprocess.run(f"samtools sort aligned_reads.sam -o Z1_sorted.bam", shell=True)
        if not r_reads:
            subprocess.run(f"{trinity_dir}/Trinity --genome_guided_bam Z1_sorted.bam --genome_guided_max_intron {max_intron} --max_memory {halfmem}G --CPU {coreamount} --seqType fq --single {f_reads} --output trinity_assembly", shell=True)
        else:
            subprocess.run(f"{trinity_dir}/Trinity --genome_guided_bam Z1_sorted.bam --genome_guided_max_intron {max_intron} --max_memory {halfmem}G --CPU 24 --seqType fq --left {f_reads} --right {r_reads} --output trinity_assembly", shell=True)
            #if len(glob.glob("trinity_assembly/Trinity-GG.fasta")) == 0:
            #    shutil.rmtree(f"{owd}/trinity_assembly")
            #    subprocess.run(f"{trinity_dir}/Trinity --genome_guided_bam Z1_sorted.bam --genome_guided_max_intron {max_intron} --max_memory {halfmem}G --CPU {coreamount} --seqType fq --left {r_reads} --right {f_reads} --output trinity_assembly", shell=True)
        with open('Time_Log.txt', 'a') as f:
            f.write(f"Trinity_Genome_Guided_Mode: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/trinity_assembly/Trinity-GG.fasta", "Assemblies_Folder/trinity_contig")
        subprocess.run(f"sed -i 's/ /_/' {owd}/Assemblies_Folder/trinity_contig", shell=True)
        shutil.rmtree(f"{owd}/trinity_assembly")
    else:
        start_time = time.time()
        try: mergedreads
        except NameError: mergedreads = None
        if mergedreads is not None:
            subprocess.run(f"{trinity_dir}/Trinity --seqType fq --max_memory {halfmem}G --single {mergedreads} --CPU {coreamount} --output trinity_assembly", shell=True)
        elif not r_reads:
            subprocess.run(f"{trinity_dir}/Trinity --seqType fq --max_memory {halfmem}G --single {f_reads} --CPU {coreamount} --output trinity_assembly", shell=True)
        else:
            subprocess.run(f"{trinity_dir}/Trinity --seqType fq --max_memory {halfmem}G --left {f_reads} --right {r_reads} --CPU {coreamount} --output trinity_assembly", shell=True)
            if len(glob.glob("trinity_assembly/Trinity.fasta")) == 0:
                shutil.rmtree(f"{owd}/trinity_assembly")
                subprocess.run(f"{trinity_dir}/Trinity --seqType fq --max_memory {halfmem}G --left {r_reads} --right {f_reads} --CPU {coreamount} --output trinity_assembly", shell=True)

        with open('Time_Log.txt', 'a') as f:
            f.write(f"Trinity_25: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/trinity_assembly/Trinity.fasta", "Assemblies_Folder/trinity_contig")
        subprocess.run(f"sed -i 's/ /_/' {owd}/Assemblies_Folder/trinity_contig", shell=True)
        shutil.rmtree(f"{owd}/trinity_assembly")

def tadpole(): #input f_reads, r_reads
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/tadpole_contig")) >= 1:
        print("Tadpole_assembly present")
    else:
        start_time = time.time()
        os.mkdir(f"{owd}/tadpole_assembly")
        os.chdir(f"{owd}/tadpole_assembly")
        if c_fig['tadpole_a'] == '1':
            adaptive_kmer()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig1.fasta k={q0} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig2.fasta k={q1} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig3.fasta k={q2} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig4.fasta k={q3} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig5.fasta k={q4} prealloc", shell=True)
            elif not r_reads:
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig1.fasta k={q0} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig2.fasta k={q1} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig3.fasta k={q2} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig4.fasta k={q3} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig5.fasta k={q4} prealloc", shell=True)
            else:
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig1.fasta k={q0} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig2.fasta k={q1} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig3.fasta k={q2} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig4.fasta k={q3} prealloc", shell=True)
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig5.fasta k={q4} prealloc", shell=True)
            subprocess.run(f"cat tadpole_contig1.fasta tadpole_contig2.fasta tadpole_contig3.fasta tadpole_contig4.fasta tadpole_contig5.fasta > tadpole_contig6.fasta", shell=True)
            subprocess.run(f"{tadpole_dir}/rename.sh in=tadpole_contig6.fasta out=tadpole_contig7.fasta prefix=CM1 ", shell=True)
            subprocess.run(f"{tadpole_dir}/sortbyname.sh in=tadpole_contig7.fasta out=tadpole_contig8.fasta length descending", shell=True)
            os.chdir(owd)
            if not os.path.exists(f"{owd}/Assemblies_Folder"):
                os.makedirs(f"{owd}/Assemblies_Folder")
            shutil.copy2(f"{owd}/tadpole_assembly/tadpole_contig8.fasta", "Assemblies_Folder/tadpole_contig")
            shutil.rmtree(f"{owd}/tadpole_assembly")

        elif c_fig['tadpole_k'] != '0':
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig1.fasta k={k0t} prealloc", shell=True)
                if k1t != 'none':
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig2.fasta k={k1t} prealloc", shell=True)
                if k2t != 'none':                    
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig3.fasta k={k2t} prealloc", shell=True)
                if k3t != 'none':                    
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig4.fasta k={k3t} prealloc", shell=True)
                if k4t != 'none':                    
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig5.fasta k={k4t} prealloc", shell=True)
            elif not r_reads:
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig1.fasta k={k0t} prealloc", shell=True)
                if k1t != 'none':
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig2.fasta k={k1t} prealloc", shell=True)
                if k2t != 'none':
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig3.fasta k={k2t} prealloc", shell=True)
                if k3t != 'none':
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig4.fasta k={k3t} prealloc", shell=True)
                if k4t != 'none':
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig5.fasta k={k4t} prealloc", shell=True)
            else:
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig1.fasta k={k0t} prealloc", shell=True)
                if k1t != 'none':
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig2.fasta k={k1t} prealloc", shell=True)
                if k2t != 'none':
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig3.fasta k={k2t} prealloc", shell=True)
                if k3t != 'none':
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig4.fasta k={k3t} prealloc", shell=True)
                if k4t != 'none':
                    subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig5.fasta k={k4t} prealloc", shell=True)
            subprocess.run(f"cat tadpole_contig* > tadpole_contig6.fasta", shell=True)
            subprocess.run(f"{tadpole_dir}/rename.sh in=tadpole_contig6.fasta out=tadpole_contig7.fasta prefix=CM1 ", shell=True)
            subprocess.run(f"{tadpole_dir}/sortbyname.sh in=tadpole_contig7.fasta out=tadpole_contig8.fasta length descending", shell=True)
            os.chdir(owd)
            if not os.path.exists(f"{owd}/Assemblies_Folder"):
                os.makedirs(f"{owd}/Assemblies_Folder")
            shutil.copy2(f"{owd}/tadpole_assembly/tadpole_contig8.fasta", "Assemblies_Folder/tadpole_contig")
            shutil.rmtree(f"{owd}/tadpole_assembly")

        else:
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={mergedreads} threads={coreamount} out=tadpole_contig1.fasta prealloc", shell=True)
            elif not r_reads:
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} threads={coreamount} out=tadpole_contig1.fasta prealloc", shell=True)
            else:
                subprocess.run(f"{tadpole_dir}/tadpole.sh in={f_reads} in2={r_reads} threads={coreamount} out=tadpole_contig1.fasta prealloc", shell=True)
            
            subprocess.run(f"{tadpole_dir}/rename.sh in=tadpole_contig1.fasta out=tadpole_contig7.fasta prefix=CM1 ", shell=True)
            subprocess.run(f"{tadpole_dir}/sortbyname.sh in=tadpole_contig7.fasta out=tadpole_contig8.fasta length descending", shell=True)
            os.chdir(owd)
            if not os.path.exists(f"{owd}/Assemblies_Folder"):
                os.makedirs(f"{owd}/Assemblies_Folder")
            shutil.copy2(f"{owd}/tadpole_assembly/tadpole_contig8.fasta", "Assemblies_Folder/tadpole_contig")
            shutil.rmtree(f"{owd}/tadpole_assembly")    

        with open('Time_Log.txt', 'a') as f:
            f.write("Tadpole: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()



def shannon(): #input f_reads, r_reads
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/shannon_contig")) >= 1:
        print("Shannon assembly already present")
    else:
        start_time = time.time()
        Temp = f"{owd}/Temp"
        try: mergedreads
        except NameError: mergedreads = None
        if mergedreads is not None:
            subprocess.run(f"{shannon_dir}/shannon_cpp shannon -o shannon_assembly -l {b1} -s {mergedreads} -k 25 -t {coreamount} -m {shannonmem} --bypass_pre_correct_read", shell=True)
        elif not r_reads:
            subprocess.run(f"{shannon_dir}/shannon_cpp shannon -o shannon_assembly -l {b1} -s {f_reads} -k 25 -t {coreamount} -m {shannonmem} --bypass_pre_correct_read", shell=True)
        else:
            subprocess.run(f"{shannon_dir}/shannon_cpp shannon -o shannon_assembly -i {b1} {b1} -p {r_reads} {f_reads} -k 25 -t {coreamount} -u 8 -m {shannonmem} --bypass_pre_correct_read", shell=True)
        with open('Time_Log.txt', 'a') as f:
            f.write("Shannon_25: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/shannon_assembly/reconstructed_seq.fasta", "Assemblies_Folder/shannon_contig")
        shutil.rmtree(f"{owd}/shannon_assembly")


def idba_tran():
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/idba_tran_contig")) >= 1:
        print("IDBA-tran assembly already present")
    else:
        start_time = time.time()  
        if c_fig['idba_tran_a'] == '1':
            adaptive_kmer()
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                shutil.copy2(f_reads, "mergedreads.gz")
                gunzip("mergedreads.gz")
                subprocess.run(f"{idba_tran_dir}/fq2fa mergedreads mergedreads.fa",shell=True)
                subprocess.run(f"{idba_tran_dir}/idba_tran -l mergedreads.fa -o idba_assembly --mink {q0} --maxk {qe} --step {qc} --num_threads {coreamount} --no_correct", shell=True)
            elif not r_reads:
                shutil.copy2(f_reads, "forward.gz")
                gunzip("forward.gz")
                subprocess.run(f"{idba_tran_dir}/fq2fa forward forward.fa",shell=True)
                subprocess.run(f"{idba_tran_dir}/idba_tran -l forward.fa -o idba_assembly --mink {q0} --maxk {qe} --step {qc} --num_threads {coreamount} --no_correct", shell=True)
            else:
                shutil.copy2(f_reads, "forward.gz")
                gunzip("forward.gz")
                shutil.copy2(r_reads, "reverse.gz")
                gunzip("reverse.gz")
                subprocess.run(f"{idba_tran_dir}/fq2fa --merge --filter forward reverse interleaved_paired.fa",shell=True)
                if b1 < 128:
                    subprocess.run(f"{idba_tran_dir}/idba_tran -r interleaved_paired.fa -o idba_assembly --mink {q0} --maxk {qe} --step {qc} --num_threads {coreamount} --no_correct", shell=True)
                else:
                    subprocess.run(f"{idba_tran_dir}/idba_tran -l interleaved_paired.fa -o idba_assembly --mink {q0} --maxk {qe} --step {qc} --num_threads {coreamount} --no_correct", shell=True)

        elif c_fig['idba_tran_k'] != '0':
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                shutil.copy2(f_reads, "mergedreads.gz")
                gunzip("mergedreads.gz")
                subprocess.run(f"{idba_tran_dir}/fq2fa mergedreads mergedreads.fa",shell=True)
                subprocess.run(f"{idba_tran_dir}/idba_tran -l mergedreads.fa -o idba_assembly --mink {k0i} --maxk {k2i} --step {k1i} --num_threads {coreamount} --no_correct", shell=True)
            elif not r_reads:
                shutil.copy2(f_reads, "forward.gz")
                gunzip("forward.gz")
                subprocess.run(f"{idba_tran_dir}/fq2fa forward forward.fa",shell=True)
                subprocess.run(f"{idba_tran_dir}/idba_tran -l forward.fa -o idba_assembly --mink {k0i} --maxk {k2i} --step {k1i} --num_threads {coreamount} --no_correct", shell=True)
            else:
                shutil.copy2(f_reads, "forward.gz")
                gunzip("forward.gz")
                shutil.copy2(r_reads, "reverse.gz")
                gunzip("reverse.gz")
                subprocess.run(f"{idba_tran_dir}/fq2fa --merge --filter forward reverse interleaved_paired.fa",shell=True)
                if b1 < 128:
                    subprocess.run(f"{idba_tran_dir}/idba_tran -r interleaved_paired.fa -o idba_assembly --mink {k0i} --maxk {k2i} --step {k1i} --num_threads {coreamount} --no_correct", shell=True)
                else:
                    subprocess.run(f"{idba_tran_dir}/idba_tran -l interleaved_paired.fa -o idba_assembly --mink {k0i} --maxk {k2i} --step {k1i} --num_threads {coreamount} --no_correct", shell=True)

        else:
            try: mergedreads
            except NameError: mergedreads = None
            if mergedreads is not None:
                shutil.copy2(f_reads, "mergedreads.gz")
                gunzip("mergedreads.gz")
                subprocess.run(f"{idba_tran_dir}/fq2fa mergedreads mergedreads.fa",shell=True)
                subprocess.run(f"{idba_tran_dir}/idba_tran -l mergedreads.fa -o idba_assembly --num_threads {coreamount} --no_correct", shell=True)
            elif not r_reads:
                shutil.copy2(f_reads, "forward.gz")
                gunzip("forward.gz")
                subprocess.run(f"{idba_tran_dir}/fq2fa forward forward.fa",shell=True)
                subprocess.run(f"{idba_tran_dir}/idba_tran -l forward.fa -o idba_assembly --num_threads {coreamount} --no_correct", shell=True)
            else:
                shutil.copy2(f_reads, "forward.gz")
                gunzip("forward.gz")
                shutil.copy2(r_reads, "reverse.gz")
                gunzip("reverse.gz")
                subprocess.run(f"{idba_tran_dir}/fq2fa --merge --filter forward reverse interleaved_paired.fa",shell=True)
                if b1 < 128:
                    subprocess.run(f"{idba_tran_dir}/idba_tran -r interleaved_paired.fa -o idba_assembly --num_threads {coreamount} --no_correct", shell=True)
                else:
                    subprocess.run(f"{idba_tran_dir}/idba_tran -l interleaved_paired.fa -o idba_assembly --num_threads {coreamount} --no_correct", shell=True)

        with open('Time_Log.txt', 'a') as f:
            f.write(f"IDBA-tran: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        subprocess.run(f"cat idba_assembly/transcript*.fa > idba_assembly/idba.fa", shell=True)
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/idba_assembly/idba.fa", "Assemblies_Folder/idba_tran_contig")
        subprocess.run(f"rm -f forward forward.fa reverse reverse.fa interleaved_paired.fa", shell=True)    
        shutil.rmtree(f"{owd}/idba_assembly")

   
def binpacker():
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/binpacker_contig")) >= 1:
        print("BinPacker assembly already present")
    else:
        start_time = time.time()
        try: mergedreads
        except NameError: mergedreads = None
        if mergedreads is not None:
            subprocess.run(f"{binpacker_dir}/BinPacker -s fq -p single -u {mergedreads} -o binpacker_assembly -k 25", shell=True)
        elif not r_reads:
            subprocess.run(f"{binpacker_dir}/BinPacker -s fq -p single -u {f_reads} -o binpacker_assembly -k 25", shell=True)
        else:
            subprocess.run(f"{binpacker_dir}/BinPacker -s fq -p pair -l {r_reads} -r {f_reads} -o binpacker_assembly -k 25", shell=True)
        with open('Time_Log.txt', 'a') as f:
            f.write("BinPacker_25: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/binpacker_assembly/BinPacker.fa", "Assemblies_Folder/binpacker_contig")
        shutil.rmtree(f"{owd}/binpacker_assembly")


def translig():
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("Assemblies_Folder/translig_contig")) >= 1:
        print("TransLig assembly already present")
    else:
        start_time = time.time()
        try: mergedreads
        except NameError: mergedreads = None
        if mergedreads is not None:
            subprocess.run(f"{translig_dir}/TransLiG -s fq -p single -u {mergedreads} -o translig_assembly -k 31", shell=True)
        elif not r_reads:
            subprocess.run(f"{translig_dir}/TransLiG -s fq -p single -u {f_reads} -o translig_assembly -k 31", shell=True)
        else:
            subprocess.run(f"{translig_dir}/TransLiG -s fq -p pair -l {r_reads} -r {f_reads} -o translig_assembly -k 31", shell=True)
        with open('Time_Log.txt', 'a') as f:
            f.write("TransLig_31: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/Assemblies_Folder"):
            os.makedirs(f"{owd}/Assemblies_Folder")
        shutil.copy2(f"{owd}/translig_assembly/TransLiG.fa", "Assemblies_Folder/translig_cont")
        os.chdir("Assemblies_Folder")
        subprocess.run(f"sed -i 's/ /_/g' translig_cont", shell=True)
        command1 = "cat translig_cont" 
        command2 = f"{seqkit_dir}/seqkit fx2tab --length" 
        command3 = """awk -F "\t" '{print $1"_length_"$4"\t"$2}'""" 
        command4 = f"{seqkit_dir}/seqkit tab2fx > translig_contig"
        subprocess.run(f"{command1}|{command2}|{command3}|{command4}", shell=True)
        subprocess.run("rm -f translig_cont", shell=True)
        os.chdir(owd)
        shutil.rmtree(f"{owd}/translig_assembly")


######################################################################################################################

######################################################################################################################
                                            ###   Post-Processing   ###
######################################################################################################################

######################################################################################################################

def transrate():
    global contig_file
    owd = os.getcwd()
    #if len(glob.glob("Assemblies_Folder/merged_contigs")) >= 1:
    #    print("TransRate merged contigs file already present")
    #else:
    #    subprocess.run(f"cat {owd}/Assemblies_Folder/*_contig > {owd}/Assemblies_Folder/merged_contigs", shell=True)        
    #merged_contigs = f"{owd}/Assemblies_Folder/merged_contigs"
    if len(glob.glob("transrate_results/merged_assemblies")) >= 1:
        print("TransRate folder already present")
        contig_file = f"{owd}/transrate_results/merged_assemblies"
    else:
        start_time = time.time()
        subprocess.run(f"sed -i 's/,//g' Assemblies_Folder/*_contig",shell=True)
        awd = f"{owd}/Assemblies_Folder"
        os.chdir(awd)
        merged_contigs = subprocess.getoutput("ls -p | grep -v / | tr '\n' ',' | sed 's/.$//' | sed 's/$/ /'")
        print(merged_contigs)        
        subprocess.run(f"{transrate_dir}/transrate --threads {coreamount} --assembly {merged_contigs} --merge-assemblies merged_assemblies", shell=True)
        contig_file = f"{owd}/transrate_results/merged_assemblies"
        shutil.move(f"{awd}/transrate_results", owd)
        os.chdir(owd)
        with open('Time_Log.txt', 'a') as f:
            f.write("TransRate_Merge: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
    os.chdir(owd)


def trim_tigs(): #input pretrim_assem = merged_assemblies, trim_len
    species_name = ntpath.basename(os.path.realpath(os.getcwd()))
    if len(glob.glob("good_merged_assemblies_busco_mod")) >= 1:
        print("Trim contigs file already present")
    else:
        start_time = time.time()
        owd = os.getcwd()
        contig_file = (os.getcwd() + "/transrate_results/merged_assemblies")
        subprocess.run(f"{tadpole_dir}/reformat.sh in={contig_file} out=good_merged_assemblies_non_sorted.fasta minlength={trim_len} fastawrap=0", shell=True)
        subprocess.run(f"{tadpole_dir}/reformat.sh in={contig_file} out=shorter_contigs.fasta maxlength={trim_len} fastawrap=0", shell=True)
        subprocess.run(f"{tadpole_dir}/sortbyname.sh in=good_merged_assemblies_non_sorted.fasta out=good_merged_assemblies.fasta length descending fastawrap=0", shell=True)
        subprocess.run(f"sed -e 's/\//_/g' -e 's/\min.*//' -e 's/\[[^][]*\]//g' good_merged_assemblies.fasta > good_merged_assemblies_busco_mod", shell=True)
        with open('Time_Log.txt', 'a') as f:
            f.write("trim_tigs: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/{species_name}_Pincho_Results_Folder"):
            os.makedirs(f"{owd}/{species_name}_Pincho_Results_Folder")
        shutil.copy2(f"good_merged_assemblies_busco_mod", f"{species_name}_Pincho_Results_Folder/{species_name}_good_merged_assemblies_busco_mod.fasta")
        shutil.copy2(f"shorter_contigs.fasta", f"{species_name}_Pincho_Results_Folder/{species_name}_shorter_contigs.fasta")

######################################################################################################################

######################################################################################################################
                                            ###   Assessment   ###
######################################################################################################################

######################################################################################################################

def hisat2():
    owd = os.getcwd()
    r_cond()
    if len(glob.glob("hisat2")) >= 1:
        print("HISAT2 folder already present.")
    else:
        start_time = time.time()
        busco_sequence = f"{owd}/good_merged_assemblies.fasta"
        os.mkdir(f"{owd}/hisat2")
        os.chdir(f"{owd}/hisat2")
        subprocess.run(f"{hisat2_dir}/hisat2-build {busco_sequence} hisat_index", shell=True)
        try: mergedreads
        except NameError: mergedreads = None
        if mergedreads is not None:
            subprocess.run(f"{hisat2_dir}/hisat2 -x hisat_index -U {mergedreads} -p {coreamount} -S aligned_reads --summary-file hisat2_summary", shell=True)
        elif not r_reads:
            subprocess.run(f"{hisat2_dir}/hisat2 -x hisat_index -U {f_reads} -p {coreamount} -S aligned_reads --summary-file hisat2_summary", shell=True)
        else:
            subprocess.run(f"{hisat2_dir}/hisat2 -x hisat_index -1 {f_reads} -2 {r_reads} -p {coreamount} -S aligned_reads --summary-file hisat2_summary --al-conc aligned_reads", shell=True)
        os.chdir(owd)
        with open('Time_Log.txt', 'a') as f:
            f.write("HISAT2: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/{species_name}_Pincho_Results_Folder"):
            os.makedirs(f"{owd}/{species_name}_Pincho_Results_Folder")
        shutil.copy2(f"hisat2/hisat2_summary", f"{species_name}_Pincho_Results_Folder/{species_name}_HISAT2_results.txt")
        shutil.rmtree(f"{owd}/hisat2")


def busco(): #input good_merged_assemblies_busco_mod
    owd = os.getcwd()
    busco_sequence = (os.getcwd() + "/good_merged_assemblies_busco_mod")
    if len(glob.glob("busco_results")) >= 1:
        print("Busco results folder already present")
    else:
        start_time = time.time()
        subprocess.run(f"python3 {busco_dir}/run_BUSCO.py -i {busco_sequence} -o busco_results -l {busco_lin} -c {coreamount} -m transcriptome --offline", shell=True)
        print(datetime.datetime.now())
        shutil.rmtree(f"{owd}/busco_downloads")
        with open('Time_Log.txt', 'a') as f:
            f.write("BUSCO: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()
        if not os.path.exists(f"{owd}/{species_name}_Pincho_Results_Folder"):
            os.makedirs(f"{owd}/{species_name}_Pincho_Results_Folder")
        os.chdir(f"{owd}/busco_results")
        subprocess.run(f"cat short_summary* > busco_results.txt", shell=True)
        os.chdir(owd)
        shutil.copy2(f"busco_results/busco_results.txt", f"{species_name}_Pincho_Results_Folder/{species_name}_busco_results.txt")
        shutil.rmtree(f"{owd}/busco_results")


def transrate_read_mapping():
    owd = os.getcwd()
    r_cond()
    assem_sequence = (os.getcwd() + "/good_merged_assemblies_busco_mod")
    if len(glob.glob("transrate_read_mapping")) >= 1:
        print("TransRate Read Mapping results folder already present")
    else:
        start_time = time.time()
        try: mergedreads
        except NameError: mergedreads = None
        if mergedreads is not None:
            subprocess.run(f"{transrate_read_mapping_dir}/transrate --assembly={assem_sequence} --left={mergedreads} --threads={coreamount} --output=transrate_read_mapping", shell=True)
        elif not r_reads:
            subprocess.run(f"{transrate_read_mapping_dir}/transrate --assembly={assem_sequence} --left={f_reads} --threads={coreamount} --output=transrate_read_mapping", shell=True)
        else:
            subprocess.run(f"{transrate_read_mapping_dir}/transrate --assembly={assem_sequence} --left={f_reads} --right={r_reads} --threads={coreamount} --output=transrate_read_mapping", shell=True)
        with open('Time_Log.txt', 'a') as f:
            f.write("TransRate_Read_Mapping: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time)))) 
            f.close()
        if not os.path.exists(f"{owd}/{species_name}_Pincho_Results_Folder"):
            os.makedirs(f"{owd}/{species_name}_Pincho_Results_Folder")
        shutil.copy2(f"transrate_read_mapping/assemblies.csv", f"{species_name}_Pincho_Results_Folder/{species_name}_transrate_results.csv")
        shutil.rmtree(f"{owd}/transrate_read_mapping")


######################################################################################################################

######################################################################################################################
                                            ###    Annotation   ###
######################################################################################################################

######################################################################################################################

#run annotation with database set in config file
def cd_hit(): #input good_merged_assemblies.fasta
    global contig_file
    start_time = time.time()
    if len(glob.glob(f"cd_contig_file")) >= 1:
        print("CD_hit already completed, skipping CD_hit")
    else:
        contig_file = (os.getcwd() + "/good_merged_assemblies.fasta")
        subprocess.run(f"{cd_hit_dir}/cd-hit-est -i {contig_file} -o cd_contig_file -c 1.00 -n 10 -M 25000 -T 0", shell=True)
        with open('Time_Log.txt', 'a') as f:
            f.write("CD_hit: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
            f.close()

def blast1():
    global contig_file, working_custom_blast_dir, working_blast_name, working_evalue, working_blast_process
    print("\n\n\nBLAST MODE 1\n\n\n")
    working_custom_blast_dir = custom_blast_dir
    working_blast_name = blast_name
    working_evalue = evalue
    working_blast_process = blast_process1
    if query_seq != '0':
        contig_file = query_seq
    elif len(glob.glob("cd_contig_file")) >= 1:
        contig_file = (os.getcwd() + "/cd_contig_file")
    else:
        contig_file = (os.getcwd() + "/good_merged_assemblies.fasta")
    blast()

def blast2():
    global contig_file, working_custom_blast_dir, working_blast_name, working_evalue, working_blast_process
    print("\n\n\nBLAST MODE 2\n\n\n")
    working_custom_blast_dir = custom_blast_dir_two
    working_blast_name = blast_name_two
    working_evalue = evalue2
    working_blast_process = blast_process2
    if query_seq != '0':
        contig_file = query_seq
    elif len(glob.glob("cd_contig_file")) >= 1:
        contig_file = (os.getcwd() + "/cd_contig_file")
    else:
        contig_file = (os.getcwd() + "/good_merged_assemblies.fasta") 
    blast()

def blast3():
    global contig_file, working_custom_blast_dir, working_blast_name, working_evalue, working_blast_process
    print("\n\n\nBLAST MODE 3: Cross Blast\n\n\n")
    contig_file = (os.getcwd() + f"/{blast_file2}_blast_results.fasta")
    working_custom_blast_dir = custom_blast_dir_three
    working_blast_name = blast_name_three
    working_evalue = evalue3
    working_blast_process = blast_process3
    blast()

def blast():
    global contig_file, species_name, blast_file2, blast_check
    species_name = ntpath.basename(os.path.realpath(os.getcwd()))    
    start_time = time.time()
    owd = os.getcwd()

    subprocess.run(f"sed -i 's/ /_/g' {contig_file}", shell=True)

    if not contig_file:
        print("No contig file was found.")
        return #stops program if contig_file is empty/not set
    else:
        print(f"\nAssembly folder found: running BLAST{working_blast_process} annotation on {contig_file} with {working_blast_name} library")

        wd = os.getcwd()
        os.chdir(working_custom_blast_dir)
        filename = ntpath.basename(f"{contig_file}")
        subprocess.run(f"{blast_dir}/bin/blast{working_blast_process} -query {contig_file} -db {working_blast_name} -evalue {working_evalue} -max_target_seqs 1 -outfmt '7 std qseqid stitle sscinames staxids' -out [{species_name}]_[{working_blast_name}]_blast_table -num_threads {coreamount}", shell=True)
        blast_file = f"[{species_name}]_[{working_blast_name}]_blast_table"
        shutil.move(blast_file, wd)
        os.chdir(wd)
        blast_file2 = f"[{species_name}]_[{working_blast_name}]_blast_table"
        print(datetime.datetime.now())
        blast_table_to_fasta()

        if not os.path.exists(f"{owd}/{species_name}_Pincho_Results_Folder"):
            os.makedirs(f"{owd}/{species_name}_Pincho_Results_Folder")
        shutil.copy2(f"{blast_file2}_blast_results.fasta", f"{species_name}_Pincho_Results_Folder/{species_name}_[{working_blast_name}]_blast_results.fasta")

    with open('Time_Log.txt', 'a') as f:
        f.write("Blastx: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
        f.close()

def blast_table_to_fasta(): #input blast table, good_merged_assemblies.fasta  
    start_time = time.time()
    subprocess.run(f"sed '/^#/d' {blast_file2} > t111", shell=True)
    subprocess.run(f"sed -i 's/^/>/' t111", shell=True)
    subprocess.run(f"cut -f1-4,11,14 t111 > t222", shell=True)
    subprocess.run(f"sed 'N;s/\\n/\\t/' {contig_file} > as111", shell=True)    
    subprocess.run(f"sed -i 's/ /_/g' t222", shell=True)
    command = """awk 'FNR==NR{a[$1]=$2;b[$1]=$3;c[$1]=$4;d[$1]=$5;e[$1]=$6;next} ($1 in a) {print $1,a[$1],b[$1],c[$1],d[$1],e[$1],$2}' t222 as111 > as222"""
    subprocess.run(command, shell=True)
    subprocess.run(f"sed 's/ /__/; s/ /__/; s/ /__/; s/ /__/; s/ /__/' as222 > as333", shell=True)
    subprocess.run(f"sed 's/ /\\n/g' as333 > {blast_file2}_blast_results.fasta", shell=True)
    subprocess.run(f"rm -f t111 t222 as111 as222 as333", shell=True) 
    with open('Time_Log.txt', 'a') as f:
        f.write("blast_conversion: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
        f.close()


######################################################################################################################

######################################################################################################################
                                            ###    Expression Analysis   ###
######################################################################################################################

######################################################################################################################


def rsem1(): #input matches.fasta, f_reads, r_reads
    global rsem_file
    r_cond()
    owd = os.getcwd()
    start_time = time.time()

    try: annotat_fasta1
    except NameError: annotat_fasta1 = None
    if annotat_fasta1 != None:
        annotation_exp = annotat_fasta1
    else:
        annotation_exp = f"{owd}/{blast_file2}_blast_results.fasta"
    if len(glob.glob("rsem1")) >= 1:
        print("RSEM1 results folder already present")
    elif not r_reads:
        os.mkdir(f"{owd}/rsem1")
        os.chdir(f"{owd}/rsem1")    
        shutil.copy2(annotation_exp, f"{owd}/rsem1/matches.fasta")
        subprocess.run(f"{trinity_dir}/util/align_and_estimate_abundance.pl --transcripts matches.fasta --seqType fq --thread_count {coreamount} --single {f_reads} --est_method RSEM --aln_method {bowtie_bin} --prep_reference --output_dir rsem_results1", shell=True)
        rsem_file = (os.getcwd() + "/rsem_results1/RSEM.genes.results")
        rsem_filter()
        os.chdir(owd)
    else:
        os.mkdir(f"{owd}/rsem1")
        os.chdir(f"{owd}/rsem1")    
        shutil.copy2(annotation_exp, f"{owd}/rsem1/matches.fasta")
        subprocess.run(f"{trinity_dir}/util/align_and_estimate_abundance.pl --transcripts matches.fasta --seqType fq --thread_count {coreamount} --left {f_reads} --right {r_reads} --est_method RSEM --aln_method {bowtie_bin} --prep_reference --output_dir rsem_results1", shell=True)
        rsem_file = (os.getcwd() + "/rsem_results1/RSEM.genes.results")
        rsem_filter()
        os.chdir(owd)
    if not os.path.exists(f"{owd}/{species_name}_Pincho_Results_Folder"):
        os.makedirs(f"{owd}/{species_name}_Pincho_Results_Folder")
    if len(glob.glob(f"{owd}/rsem1/rsem_results1/RSEM.genes.results")) >= 1:
        shutil.copy2(f"{owd}/rsem1/rsem_results1/RSEM.genes.results", f"{species_name}_Pincho_Results_Folder/{species_name}_RSEM1.genes.results")
    with open('Time_Log.txt', 'a') as f:
        f.write("RSEM1: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
        f.close()

def rsem2(): #input matches.fasta, f_reads, r_reads
    global rsem_file
    r_cond()
    owd = os.getcwd()
    start_time = time.time()

    try: annotat_fasta2
    except NameError: annotat_fasta2 = None
    if annotat_fasta2 != None:
        annotation_exp = annotat_fasta2
    else:
        annotation_exp = f"{owd}/{blast_file2}_blast_results.fasta"
    if len(glob.glob("rsem2")) >= 1:
        print("RSEM2 results folder already present")
    elif not r_reads:
        os.mkdir(f"{owd}/rsem2")
        os.chdir(f"{owd}/rsem2")    
        shutil.copy2(annotation_exp, f"{owd}/rsem2/matches.fasta")
        subprocess.run(f"{trinity_dir}/util/align_and_estimate_abundance.pl --transcripts matches.fasta --seqType fq --thread_count {coreamount} --single {f_reads} --est_method RSEM --aln_method {bowtie_bin} --prep_reference --output_dir rsem_results2", shell=True)
        rsem_file = (os.getcwd() + "/rsem_results2/RSEM.genes.results")
        rsem_filter()
        os.chdir(owd)
    else:
        os.mkdir(f"{owd}/rsem2")
        os.chdir(f"{owd}/rsem2")    
        shutil.copy2(annotation_exp, f"{owd}/rsem2/matches.fasta")
        subprocess.run(f"{trinity_dir}/util/align_and_estimate_abundance.pl --transcripts matches.fasta --seqType fq --thread_count {coreamount} --left {f_reads} --right {r_reads} --est_method RSEM --aln_method {bowtie_bin} --prep_reference --output_dir rsem_results2", shell=True)
        rsem_file = (os.getcwd() + "/rsem_results2/RSEM.genes.results")
        rsem_filter()
        os.chdir(owd)
    if not os.path.exists(f"{owd}/{species_name}_Pincho_Results_Folder"):
        os.makedirs(f"{owd}/{species_name}_Pincho_Results_Folder")
    if len(glob.glob(f"{owd}/rsem2/rsem_results2/RSEM.genes.results")) >= 1:
        shutil.copy2(f"{owd}/rsem2/rsem_results2/RSEM.genes.results", f"{species_name}_Pincho_Results_Folder/{species_name}_RSEM2.genes.results")
    with open('Time_Log.txt', 'a') as f:
        f.write("RSEM2: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
        f.close()

def rsem_filter(): #input rsem_file
    subprocess.run(f"sed '1d' {rsem_file} > t111", shell=True) #remove row 1
    subprocess.run(f"cut -f1,3-7 t111 > t222", shell=True) #cut and keep col 1, 3-7
    subprocess.run(f"sed -i 's/__/\t/g' t222", shell=True) #replace double underscore with tab
    subprocess.run(f"cut -f1-2,6-11 t222 > t333", shell=True)
    command = """awk 'BEGIN{FS=OFS="\t"} {$2 = $2 OFS $2} 1' t333 > t444"""
    subprocess.run(command, shell=True)
    command = """awk 'cnt[$3]++{$3=$3"."cnt[$3]-1} 1' t444 > t555"""
    subprocess.run(command, shell=True)
    subprocess.run(f"sed -i 's/ /\t/g' t555", shell=True)
    command = """awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $4}' t555 > RSEM_FINAL"""
    subprocess.run(command, shell=True)
    subprocess.run(f"""sed -i 1i"gene_id\tuniprot_id\ttranscript_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tblast_match" RSEM_FINAL_RESULTS""", shell=True)
    subprocess.run(f"rm -f t111 t222 t333 t444 t555", shell=True)
    print("RSEM data filtered as RSEM_FINAL_RESULTS")


def kallisto1(): #input matches.fasta, f_reads, r_reads
    start_time = time.time()
    r_cond()
    owd = os.getcwd()
    try: annotat_fasta1
    except NameError: annotat_fasta1 = None
    if annotat_fasta1 != None:
        annotation_exp = annotat_fasta1
    else:
        annotation_exp = f"{owd}/{blast_file2}_blast_results.fasta"
    if len(glob.glob("kallisto1")) >= 1:
        print("Kallisto1 results folder already present")
    elif not r_reads:
        os.mkdir(f"{owd}/kallisto1")
        os.chdir(f"{owd}/kallisto1")
        subprocess.run(f"{kallisto_dir}/kallisto index -i transcripts.idx {annotation_exp}", shell=True)
        subprocess.run(f"{kallisto_dir}/kallisto quant -i transcripts.idx -o output1 -b 100 --single -l 100 -s 20 {f_reads} -t {coreamount}", shell=True)
        os.chdir(owd)
    else:
        os.mkdir(f"{owd}/kallisto1")
        os.chdir(f"{owd}/kallisto1")
        subprocess.run(f"{kallisto_dir}/kallisto index -i transcripts.idx {annotation_exp}", shell=True)
        subprocess.run(f"{kallisto_dir}/kallisto quant -i transcripts.idx -o output1 -b 100 {f_reads} {r_reads} -t {coreamount}", shell=True)
        os.chdir(owd)
    if not os.path.exists(f"{owd}/{species_name}_Pincho_Results_Folder"):
        os.makedirs(f"{owd}/{species_name}_Pincho_Results_Folder")
    shutil.copy2(f"{owd}/kallisto1/output1/abundance.tsv", f"{species_name}_Pincho_Results_Folder/{species_name}_kallisto1_abundance.tsv")
    with open('Time_Log.txt', 'a') as f:
        f.write("Kallisto1: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
        f.close()

def kallisto2(): #input matches.fasta, f_reads, r_reads
    start_time = time.time()
    r_cond()
    owd = os.getcwd()
    try: annotat_fasta2
    except NameError: annotat_fasta2 = None
    if annotat_fasta2 != None:
        annotation_exp = annotat_fasta2
    else:
        annotation_exp = f"{owd}/{blast_file2}_blast_results.fasta"
    if len(glob.glob("kallisto2")) >= 1:
        print("Kallisto2 results folder already present")
    elif not r_reads:
        os.mkdir(f"{owd}/kallisto2")
        os.chdir(f"{owd}/kallisto2")
        subprocess.run(f"{kallisto_dir}/kallisto index -i transcripts.idx {annotation_exp}", shell=True)
        subprocess.run(f"{kallisto_dir}/kallisto quant -i transcripts.idx -o output2 -b 100 --single -l 100 -s 20 {f_reads} -t {coreamount}", shell=True)
        os.chdir(owd)
    else:
        os.mkdir(f"{owd}/kallisto2")
        os.chdir(f"{owd}/kallisto2")
        subprocess.run(f"{kallisto_dir}/kallisto index -i transcripts.idx {annotation_exp}", shell=True)
        subprocess.run(f"{kallisto_dir}/kallisto quant -i transcripts.idx -o output2 -b 100 {f_reads} {r_reads} -t {coreamount}", shell=True)
        os.chdir(owd)
    if not os.path.exists(f"{owd}/{species_name}_Pincho_Results_Folder"):
        os.makedirs(f"{owd}/{species_name}_Pincho_Results_Folder")
    shutil.copy2(f"{owd}/kallisto2/output2/abundance.tsv", f"{species_name}_Pincho_Results_Folder/{species_name}_kallisto2_abundance.tsv")
    with open('Time_Log.txt', 'a') as f:
        f.write("Kallisto2: %s \n" % str(datetime.timedelta(seconds=round(time.time() - start_time))))
        f.close()

def cleanup():
    subprocess.run(f"rm -f forward_paired.fastq.gz forward_unpaired.fastq.gz reverse_paired.fastq.gz reverse_unpaired.fastq.gz forward_paired.cor.fq.gz reverse_paired.cor.fq.gz good_merged_assemblies_non_sorted.fasta hisat2/aligned_reads* hisat2/hisat_index* rsem/rsem_results/bowtie2* length.txt good_merged_assemblies* cd_contig_file.clstr shorter_contigs.fasta", shell=True)    
    if len(glob.glob("transrate_results")) >= 1:
        shutil.rmtree("transrate_results")


######################################################################################################################

######################################################################################################################
                                                ###   END April 2 2021
######################################################################################################################

######################################################################################################################
