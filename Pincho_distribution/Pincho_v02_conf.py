#!/usr/bin/env python3
import tkinter
import subprocess
from tkinter import *
import tkinter.filedialog
import tkinter.filedialog as filedialog
from tkinter import ttk
import os
from os import path
from multiprocessing import Process
import shutil
import sys

### master directories
cwd = os.getcwd()
bin_dir = os.path.join(sys.path[0], 'bin')
lib_dir = os.path.join(sys.path[0], 'lib')
pinc_dir = sys.path[0]
conf_file = os.path.join(sys.path[0], 'pincho_v02.conf')
dirimage = f"{pinc_dir}/directory_image.png"
ancestral_file = os.path.join(lib_dir, "eukaryota_odb10.2019-11-20/eukaryota_odb10/ancestral")
adapter_file = os.path.join(bin_dir, "Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa")
blastx_db_file = os.path.join(lib_dir, "uniprot_db/uniprot_sprot.fasta")
blastn_db_file = os.path.join(lib_dir, "kegg_libraries/xla_npr_kegg/x_n_kegg_db.nhr")
blastx_db_file2 = os.path.join(lib_dir, "trembl_amphi/trembl_amphibia.fasta")

def browsefunc1():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v13.set(filename)

def browsefunc2():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v14.set(filename)

def browsefunc3():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v18.set(filename)

def browsefunc4():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v19.set(filename)

def browsefunc5():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v22.set(filename)

def browsefunc6():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v23.set(filename)

def browsefunc7():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v24.set(filename)

def browsefunc8():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v25.set(filename)

def browsefunc9():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v27.set(filename)

def browsefunc10():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v29.set(filename)

def browsefunc11():
    filename = filedialog.askopenfilename()
    pathlabel.config(text=filename)
    v33.set(filename)

def start():
	owd = os.getcwd()
	if not v2.get():
		v2.set('"1e-10"')
	if not v3.get():
		v3.set('"1e-10"')
	if not v4.get():
		v4.set(0)
	if not v5.get():
		v5.set(0)
	if not v6.get():
		v6.set(0)
	if not v7.get():
		v7.set(0)
	if not v8.get():
		v8.set(0)
	if not v9.get():
		v9.set(0)
	if not v10.get():
		v10.set(0)
	if not v11.get():
		v11.set(0)
	if not v12.get():
		v12.set(0)
	if not v13.get():
		v13.set(0)
	if not v14.get():
		v14.set(0)
	if not v15.get():
		v15.set(0)
	if not v16.get():
		v16.set(0)
	if not v17.get():
		v17.set(0)
	if not v18.get():
		v18.set(0)
	if not v19.get():
		v19.set(0)
	if not v20.get():
		v20.set(0)
	if not v21.get():
		v21.set(0)
	if not v22.get():
		v22.set(0)
	if not v23.get():
		v23.set(0)
	if not v24.get():
		v24.set(0)
	if not v25.get():
		v25.set(0)
	if not v26.get():
		v26.set(0)
	if not v27.get():
		v27.set(0)
	if not v28.get():
		v28.set(0)
	if not v29.get():
		v29.set(0)
	if not v30.get():
		v30.set(0)
	if not v31.get():
		v31.set(0)
	if not v32.get():
		v32.set(0)
	if not v33.get():
		v33.set(0)
	if not sal1.get():
		sal1.set(0)
	if not sal2.get():
		sal2.set(0)
	if not blast_resource.get():
		blast_resource.set(0)
	if not v34.get():
		v34.set(0)
	conf = ("trimmomatic =\t%d\nrcorrector =\t%d\nabyss \t=\t%d\nspades =\t%d\nmegahit =\t%d\ntransabyss =\t%d\nrnaspades =\t%d\ntadpole =\t%d\noases \t=\t%d\nshannon =\t%d\ntrinity =\t%d\nidba_tran =\t%d\ntranslig =\t%d\nbinpacker =\t%d\ntransrate  =\t%d\ntrim_tigs =\t%d\nbusco \t=\t%d\ntransrate2 =\t%d\nhisat2 =\t%d\ncd_hit =\t%d\nblast \t=\t%d\nrsem1 \t=\t%d\nkallisto1 =\t%d\nblast2 =\t%d\nabyss_a =\t%d\nspades_a =\t%d\nmegahit_a =\t%d\ntransabyss_a =\t%d\nrnaspades_a =\t%d\ntadpole_a =\t%d\noases_a =\t%d\nidba_tran_a =\t%d\ntrim_length =\t%d\nblast_evalue1 =\t%s\nblast_evalue3 =\t%s\nabyss_k =\t%s\nspades_k =\t%s\nmegahit_k =\t%s\ntransabyss_k =\t%s\nrnaspades_k =\t%s\ntadpole_k =\t%s\noases_k =\t%s\nidba_tran_k =\t%s\nSRR \t=\t%s\nread1 \t=\t%s\nread2 \t=\t%s\nnfolders =\t%d\nnthreads =\t%d\nnmemory =\t%d\ncleanup \t=\t%d\nadapter_seqs \t=\t%s\nbusco_lin \t=\t%s\nblast_type1 \t=\t%s\nblast_type3 \t=\t%s\nblast_db1 \t=\t%s\nblast_db3 \t=\t%s\nquery_seq \t=\t%s\nannotation_fasta1 \t=\t%s\nblast_type2 \t=\t%s\nblast_db2 \t=\t%s\nblast_evalue2 \t=\t%s\ngenome_guide \t=\t%s\nmax_intron \t=\t%s\nrsem2 \t=\t%s\nkallisto2 \t=\t%s\nannotation_fasta2 \t=\t%s\nsalmon1 \t=\t%s\nsalmon2 \t=\t%s\nnthreads_blast =\t%d\nSRR_folder_names \t=\t%s" % \
		(CV1.get(), CV2.get(), CV3.get(), CV4.get(), CV5.get(), CV6.get(), CV7.get(), CV8.get(), CV9.get(), CV10.get(), CV11.get(), CV12.get(), CV13.get(), CV14.get(), CV16.get(), CV17.get(), CV18.get(), CV19.get(), CV20.get(), CV22.get(), CV23.get(), CV24.get(), CV25.get(), CV27.get(), CV3b.get(), CV4b.get(), CV5b.get(), CV6b.get(), CV7b.get(), CV8b.get(), CV9b.get(), CV12b.get(), v1.get(), v2.get(), v3.get(), v4.get(), v5.get(), v6.get(), v7.get(), v8.get(), v9.get(), v10.get(), v11.get(), v12.get(), v13.get(), v14.get(), v15.get(), v16.get(), v17.get(), CV28.get(), v18.get(), v19.get(), v20.get(), v21.get(), v22.get(), v23.get(), v24.get(), v25.get(), v26.get(), v27.get(), v28.get(), v29.get(), v30.get(), v31.get(), v32.get(), v33.get(), sal1.get(), sal2.get(), blast_resource.get(), v34.get()))
	os.chdir(pinc_dir)
	with open('pincho_v02.conf', 'w+') as f:
		f.write(conf)
		f.close()
	os.chdir(owd)
	top.wm_state('iconic')
	run_pincho()
	quit()

def run_pincho():
	subprocess.run(f"python3 {pinc_dir}/Pincho_v02.py -conf {pinc_dir}/pincho_v02.conf", shell=True)


def high_quality():
	CV1.set(1)
	CV2.set(1)
	CV3.set(0)
	CV4.set(0)
	CV5.set(0)
	CV6.set(1)
	CV7.set(1)
	CV8.set(0)
	CV9.set(0)
	CV10.set(0)
	CV11.set(0)
	CV12.set(0)
	CV13.set(1)
	CV14.set(0)
	CV3b.set(0)
	CV4b.set(0)
	CV5b.set(0)
	CV6b.set(1)
	CV7b.set(1)
	CV8b.set(0)
	CV9b.set(0)
	CV10b.set(0)
	CV11b.set(0)
	CV12b.set(0)
	CV13b.set(1)
	CV14b.set(0)
	CV15.set(1)
	CV16.set(1)
	CV17.set(1)
	CV18.set(1)
	CV19.set(0)
	CV20.set(0)
	CV22.set(1)
	CV23.set(1)
	CV24.set(1)
	CV25.set(1)
	CV27.set(1)
	CV28.set(1)
	v1.set(0)
	v2.set('"1e-10"')
	v3.set('"1e-10"')
	v4.set("")
	v5.set("")
	v6.set("")
	v7.set("")
	v8.set("")
	v9.set("")
	v10.set("")
	v11.set("")
	v12.set("")
	v13.set("")
	v14.set("")
	v15.set(1)
	v16.set(0)
	v17.set(180)
	v18.set(adapter_file)
	v19.set(ancestral_file)
	v20.set("mmseqs")
	v21.set("n")
	v22.set(blastx_db_file2)
	v23.set(blastn_db_file)
	v24.set("")
	v25.set("")
	v26.set("mmseqs")
	v27.set(blastx_db_file)
	v28.set('"1e-10"')
	v29.set("")
	v30.set(0)
	v31.set(1)
	v32.set(1)
	v33.set("")
	sal1.set(1)
	sal2.set(1)
	blast_resource.set(0)
	v34.set("")	

def clear():
	CV1.set(0)
	CV2.set(0)
	CV3.set(0)
	CV4.set(0)
	CV5.set(0)
	CV6.set(0)
	CV7.set(0)
	CV8.set(0)
	CV9.set(0)
	CV10.set(0)
	CV11.set(0)
	CV12.set(0)
	CV13.set(0)
	CV14.set(0)
	CV3b.set(0)
	CV4b.set(0)
	CV5b.set(0)
	CV6b.set(0)
	CV7b.set(0)
	CV8b.set(0)
	CV9b.set(0)
	CV10b.set(0)
	CV11b.set(0)
	CV12b.set(0)
	CV13b.set(0)
	CV14b.set(0)
	CV15.set(0)
	CV16.set(0)
	CV17.set(0)
	CV18.set(0)
	CV19.set(0)
	CV20.set(0)
	CV22.set(0)
	CV23.set(0)
	CV24.set(0)
	CV25.set(0)
	CV27.set(0)
	CV28.set(0)
	v1.set(0)
	v2.set('"1e-10"')
	v3.set('"1e-10"')
	v4.set("")
	v5.set("")
	v6.set("")
	v7.set("")
	v8.set("")
	v9.set("")
	v10.set("")
	v11.set("")
	v12.set("")
	v13.set("")
	v14.set("")
	v15.set(1)
	v16.set(0)
	v17.set(4)
	v18.set("")
	v19.set("")
	v20.set("")
	v21.set("")
	v22.set("")
	v23.set("")
	v24.set("")
	v25.set("")
	v26.set("")
	v27.set("")
	v28.set('"1e-10"')
	v29.set("")
	v30.set(0)
	v31.set(0)
	v32.set(0)
	v33.set("")	
	sal1.set(0)
	sal2.set(0)
	blast_resource.set(0)
	v34.set("")	

def test():
	CV1.set(1)
	CV2.set(0)
	CV3.set(0)
	CV4.set(0)
	CV5.set(1)
	CV6.set(0)
	CV7.set(0)
	CV8.set(0)
	CV9.set(0)
	CV10.set(0)
	CV11.set(0)
	CV12.set(0)
	CV13.set(0)
	CV14.set(0)
	CV3b.set(1)
	CV4b.set(1)
	CV5b.set(1)
	CV6b.set(1)
	CV7b.set(1)
	CV8b.set(1)
	CV9b.set(1)
	CV10b.set(1)
	CV11b.set(1)
	CV12b.set(1)
	CV13b.set(1)
	CV14b.set(1)
	CV15.set(0)
	CV16.set(1)
	CV17.set(1)
	CV18.set(0)
	CV19.set(0)
	CV20.set(0)
	CV22.set(0)
	CV23.set(1)
	CV24.set(1)
	CV25.set(1)
	CV27.set(1)
	CV28.set(1)
	v1.set(0)
	v2.set('"1e-10"')
	v3.set('"1e-10"')
	v4.set("")
	v5.set("")
	v6.set("")
	v7.set("")
	v8.set("")
	v9.set("")
	v10.set("")
	v11.set("")
	v12.set("SRR11880886")
	v13.set("")
	v14.set("")
	v15.set(1)
	v16.set(0)
	v17.set(180)
	v18.set(adapter_file)
	v19.set(ancestral_file)
	v20.set("mmseqs")
	v21.set("n")
	v22.set(blastx_db_file2)
	v23.set(blastn_db_file)
	v24.set("")
	v25.set("")
	v26.set("mmseqs")
	v27.set(blastx_db_file)
	v28.set('"1e-10"')
	v29.set("")
	v30.set(0)
	v31.set(1)
	v32.set(1)
	v33.set("")	
	sal1.set(1)
	sal2.set(1)
	blast_resource.set(0)
	v34.set("SRR11880886_Test_Species")	

top = tkinter.Tk()
top.iconphoto(False, tkinter.PhotoImage(file=f'{pinc_dir}/pincho_icon.png'))
top.title("Pincho v0.1")
top.call('tk', 'scaling', 0.1)
CV1 = IntVar()
CV2 = IntVar()
CV3 = IntVar()
CV4 = IntVar()
CV5 = IntVar()
CV6 = IntVar()
CV7 = IntVar()
CV8 = IntVar()
CV9 = IntVar()
CV10 = IntVar()
CV11 = IntVar()
CV12 = IntVar()
CV13 = IntVar()
CV14 = IntVar()
CV3b = IntVar()
CV4b = IntVar()
CV5b = IntVar()
CV6b = IntVar()
CV7b = IntVar()
CV8b = IntVar()
CV9b = IntVar()
CV10b = IntVar()
CV11b = IntVar()
CV12b = IntVar()
CV13b = IntVar()
CV14b = IntVar()
CV15 = IntVar()
CV16 = IntVar()
CV17 = IntVar()
CV18 = IntVar()
CV19 = IntVar()
CV20 = IntVar()
CV22 = IntVar()
CV23 = IntVar()
CV24 = IntVar()
CV25 = IntVar()
CV27 = IntVar()
CV28 = IntVar()
sal1 = IntVar()
sal2 = IntVar()
blast_resource = IntVar()

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
v12 = StringVar()
tkinter.Label(top, text = "Single Run:", fg="blue").grid(row = 70, column = 17, sticky=W) 
tkinter.Label(top, text = "SRA Accession #").grid(row = 80, column = 17) 
tkinter.Entry(top, text=v12).grid(row = 80, column = 18) 
v12.set("")

v34 = StringVar()
tkinter.Label(top, text = "If SRA, create folder names:").grid(row = 81, column = 17) 
tkinter.Entry(top, text=v34).grid(row = 81, column = 18) 
v34.set("")

tkinter.Label(top, text = "OR").grid(row = 90, column = 17) 

v13 = StringVar()
tkinter.Label(top, text = "Read 1 (Forward)").grid(row = 100, column = 17) 
tkinter.Entry(top, text=v13).grid(row = 100, column = 18) 
v13.set("")
browsebutton = Button(top, text="Browse", command=browsefunc1).grid(row = 100, column = 19, sticky=W)
pathlabel = Label(top)

v14 = StringVar()
tkinter.Label(top, text = "Read 2 (Reverse)").grid(row = 110, column = 17) 
tkinter.Entry(top, text=v14).grid(row = 110, column = 18) 
v14.set("")
browsebutton = Button(top, text="Browse", command=browsefunc2).grid(row = 110, column = 19, sticky=W)
pathlabel = Label(top)
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
v15 = IntVar()
tkinter.Label(top, text = "Bulk Run:", fg="blue").grid(row = 120, column = 17, sticky=W) 
tkinter.Label(top, text = "# of Folders", fg="red").grid(row = 130, column = 17) 
tkinter.Entry(top, text=v15).grid(row = 130, column = 18) 
v15.set(1)

tkinter.Label(top, text = "Directory must follow the structure below:", fg="green").grid(row = 140, column = 17, columnspan =2) 
tkinter.Label(top, text = "Working folders must be inside same parent", fg="green").grid(row = 150, column = 17, columnspan =2) 
tkinter.Label(top, text = "No Folder Limit. All folders will be processed.", fg="green").grid(row = 160, column = 17, columnspan =2) 
photo = PhotoImage(file = dirimage)
tkinter.Label(top, image = photo).grid(row = 170, rowspan =100, column = 17, columnspan =2) 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
tkinter.Label(top, text = "Pre-Processing:", fg="blue").grid(row = 0, column = 0, sticky=W)

tkinter.Checkbutton(top, text = "Trimmomatic", variable = CV1, onvalue = 1, offvalue=0).grid(row = 10,sticky=W)

v18 = StringVar()
tkinter.Label(top, text = "Adapter Sequence & Poly-A Removal:").grid(row = 10, column = 1) 
tkinter.Entry(top, text=v18).grid(row = 10, column = 2) 
v18.set("adapter fasta")
browsebutton = Button(top, text="Browse", command=browsefunc3).grid(row = 10, column = 3, sticky=W)
pathlabel = Label(top)

tkinter.Checkbutton(top, text = "Rcorrector", variable = CV2, onvalue = 1, offvalue =0).grid(row = 20,sticky=W)
tkinter.Label(top, text = "Error Correction for Illumina NGS Data").grid(row = 20, column = 1) 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
tkinter.Label(top, text = "Assembly:", fg="blue").grid(row = 30, column = 0, sticky=W)

tkinter.Checkbutton(top, text = "ABySS", variable = CV3, onvalue = 1, offvalue=0).grid(row = 50,sticky=W)
tkinter.Checkbutton(top, text = "SPAdes", variable = CV4, onvalue = 1, offvalue=0).grid(row = 60,sticky=W)
tkinter.Checkbutton(top, text = "MEGAHIT", variable = CV5, onvalue = 1, offvalue=0).grid(row = 70,sticky=W)
tkinter.Checkbutton(top, text = "trans-ABySS", variable = CV6, onvalue = 1, offvalue=0).grid(row = 80,sticky=W)
tkinter.Checkbutton(top, text = "rnaSPAdes", variable = CV7, onvalue = 1, offvalue=0).grid(row = 90,sticky=W)
tkinter.Checkbutton(top, text = "Tadpole", variable = CV8, onvalue = 1, offvalue=0).grid(row = 100,sticky=W)
tkinter.Checkbutton(top, text = "Oases", variable = CV9, onvalue = 1, offvalue=0).grid(row = 170,sticky=W)
tkinter.Checkbutton(top, text = "Shannon Cpp", variable = CV10, onvalue = 1, offvalue=0).grid(row = 120,sticky=W)
tkinter.Checkbutton(top, text = "Trinity", variable = CV11, onvalue = 1, offvalue=0).grid(row = 130,sticky=W)
tkinter.Checkbutton(top, text = "IDBA-tran", variable = CV12, onvalue = 1, offvalue=0).grid(row = 180,sticky=W)
tkinter.Checkbutton(top, text = "TransLig", variable = CV13, onvalue = 1, offvalue=0).grid(row = 160,sticky=W)
tkinter.Checkbutton(top, text = "BinPacker", variable = CV14, onvalue = 1, offvalue=0).grid(row = 110,sticky=W)

tkinter.Label(top, text = "k-list (max:5) k-step (max:3, must be odd)").grid(row = 40, column = 1) 
tkinter.Label(top, text = "adaptive k-mers").grid(row = 40, column = 2) 


v4 = StringVar()
tkinter.Entry(top, text=v4).grid(row = 50, column = 1) 
v4.set("21,41,61,81,99")
tkinter.Checkbutton(top, text = "", variable = CV3b, onvalue = 1, offvalue =0).grid(row = 50,column = 2,sticky=W)

v5 = StringVar()
tkinter.Entry(top, text=v5).grid(row = 60, column = 1) 
v5.set("21,41,61,81,99")
tkinter.Checkbutton(top, text = "", variable = CV4b, onvalue = 1, offvalue =0).grid(row = 60,column = 2,sticky=W)

v6 = StringVar()
tkinter.Entry(top, text=v6).grid(row = 70, column = 1) 
tkinter.Checkbutton(top, text = "", variable = CV5b, onvalue = 1, offvalue =0).grid(row = 70, column = 2,sticky=W)
v6.set("21,41,61,81,99")

v7 = StringVar()
tkinter.Entry(top, text=v7).grid(row = 80, column = 1) 
v7.set("21,41,61,81,99")
tkinter.Checkbutton(top, text = "", variable = CV6b, onvalue = 1, offvalue=0).grid(row = 80, column = 2,sticky=W)

v8 = StringVar()
tkinter.Entry(top, text=v8).grid(row = 90, column = 1) 
v8.set("21,41,61,81,99")
tkinter.Checkbutton(top, text = "", variable = CV7b, onvalue = 1, offvalue=0).grid(row = 90, column = 2,sticky=W)

v9 = StringVar()
tkinter.Entry(top, text=v9).grid(row = 100, column = 1) 
v9.set("21,41,61,81,99")
tkinter.Checkbutton(top, text = "", variable = CV8b, onvalue = 1, offvalue =0).grid(row = 100, column = 2,sticky=W)

tkinter.Label(top, text = "25").grid(row = 110, column = 1) 

tkinter.Label(top, text = "25").grid(row = 120, column = 1)

tkinter.Label(top, text = "25").grid(row = 130, column = 1) 
tkinter.Label(top, text = "Genome Guided Mode Available!   >>>>>", fg="gray").grid(row = 140, column = 1) 
tkinter.Label(top, text = "Genome max intron Length            >>>>>", fg="gray").grid(row = 150, column = 1) 
v30 = StringVar()
tkinter.Entry(top, text=v30).grid(row = 150, column = 2) 
v30.set(10000)
tkinter.Label(top, text = "bp").grid(row = 150, column = 3, sticky=W) 

v29 = StringVar()
tkinter.Entry(top, text=v29).grid(row = 140, column = 2) 
v29.set("genome fasta")
browsebutton = Button(top, text="Browse", command=browsefunc10).grid(row = 140, column = 3, sticky=W)
pathlabel = Label(top)


tkinter.Label(top, text = "31").grid(row = 160, column = 1) 

v10 = StringVar()
tkinter.Entry(top, text=v10).grid(row = 170, column = 1) 
v10.set("21,10,51")
tkinter.Checkbutton(top, text = "", variable = CV9b, onvalue = 1, offvalue =0).grid(row = 170, column = 2,sticky=W)

v11 = StringVar()
tkinter.Entry(top, text=v11).grid(row = 180, column = 1) 
v11.set("21,10,51")
tkinter.Checkbutton(top, text = "", variable = CV12b, onvalue = 1, offvalue =0).grid(row = 180, column = 2,sticky=W)
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
tkinter.Label(top, text = "Post-Processing:", fg="blue").grid(row = 190, column = 0, sticky=W)

tkinter.Checkbutton(top, text = "TransRate", variable = CV16, onvalue = 1, offvalue=0).grid(row = 200,sticky=W)
tkinter.Label(top, text = "Prepares Assembly for Annotation").grid(row = 200, column = 1) 

tkinter.Checkbutton(top, text = "Trim Contigs", variable = CV17, onvalue = 1, offvalue=0).grid(row = 210,sticky=W)

v1 = IntVar()
tkinter.Label(top, text = "Remove transcripts under:").grid(row = 210, column = 1)
tkinter.Entry(top, text=v1).grid(row = 210, column = 2) 
tkinter.Label(top, text = "bp").grid(row = 210, column = 3, sticky=W) 
v1.set(300)

tkinter.Checkbutton(top, text = "CD-HIT", variable = CV22, onvalue = 1, offvalue=0).grid(row = 220,sticky=W)
tkinter.Label(top, text = "Remove Duplicate Sequences before Annotation").grid(row = 220, column = 1) 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
tkinter.Label(top, text = "Assembly Assessment:", fg="blue").grid(row = 230, column = 0, sticky=W)

tkinter.Checkbutton(top, text = "BUSCO", variable = CV18, onvalue = 1, offvalue=0).grid(row = 240,sticky=W)

v19 = StringVar()
tkinter.Label(top, text = "BUSCO Lineage File:").grid(row = 240, column = 1) 
tkinter.Entry(top, text=v19).grid(row = 240, column = 2) 
v19.set("busco ancestral file")
browsebutton = Button(top, text="Browse", command=browsefunc4).grid(row = 240, column = 3, sticky=W)
pathlabel = Label(top)

tkinter.Checkbutton(top, text = "TransRate Assessment", variable = CV19, onvalue = 1, offvalue=0).grid(row = 250,sticky=W)
tkinter.Label(top, text = "TransRate Score i.e. N50/N90").grid(row = 250, column = 1) 

tkinter.Checkbutton(top, text = "HISAT2", variable = CV20, onvalue = 1, offvalue=0).grid(row = 260,sticky=W)
tkinter.Label(top, text = "Raw Data Utilization %").grid(row = 260, column = 1) 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
tkinter.Label(top, text = "Annotation:", fg="blue").grid(row = 270, column = 0, sticky=W)
tkinter.Label(top, text = "If Starting with Annotation Provide Fasta to Annotate:", fg="gray").grid(row = 270, column = 1)
v24 = StringVar()
tkinter.Entry(top, text=v24).grid(row = 270, column = 2) 
v24.set("")
browsebutton = Button(top, text="Browse", command=browsefunc7).grid(row = 270, column = 3, sticky=W)
pathlabel = Label(top)

#tkinter.Checkbutton(top, text = "BLAST", variable = CV23, onvalue = 1, offvalue=0).grid(row = 260,sticky=W)
tkinter.Label(top, text = "Blast1: type x or n or p or mmseqs").grid(row = 280, column = 0, sticky=W)

v20 = StringVar()
tkinter.Entry(top, text=v20).grid(row = 280, column = 1) 
v20.set("")

v22 = StringVar()
tkinter.Label(top, text = "db file (or fasta to create new db)      >>>>>").grid(row = 290, column = 1) 
tkinter.Entry(top, text=v22).grid(row = 290, column = 2) 
v22.set("db file (or fasta)")
browsebutton = Button(top, text="Browse", command=browsefunc5).grid(row = 290, column = 3, sticky=W)
pathlabel = Label(top)


tkinter.Label(top, text = "Blast2: type x or n or p or mmseqs").grid(row = 300, column = 0, sticky=W)

v26 = StringVar()
tkinter.Entry(top, text=v26).grid(row = 300, column = 1) 
v26.set("")

v27 = StringVar()
tkinter.Label(top, text = "db file (or fasta to create new db)      >>>>>").grid(row = 310, column = 1) 
tkinter.Entry(top, text=v27).grid(row = 310, column = 2) 
v27.set("db file (or fasta)")
browsebutton = Button(top, text="Browse", command=browsefunc9).grid(row = 310, column = 3, sticky=W)
pathlabel = Label(top)


#tkinter.Checkbutton(top, text = "Cross_BLAST", variable = CV27, onvalue = 1, offvalue=0).grid(row = 280,sticky=W)
tkinter.Label(top, text = "Cross Blast2: type x or n or p or mmseqs").grid(row = 320, column = 0, sticky=W)

v21 = StringVar()
tkinter.Entry(top, text=v21).grid(row = 320, column = 1) 
v21.set("")

v23 = StringVar()
tkinter.Label(top, text = "db file (or fasta to create new db)      >>>>>").grid(row = 330, column = 1) 
tkinter.Entry(top, text=v23).grid(row = 330, column = 2) 
v23.set("db file (or fasta)")
browsebutton = Button(top, text="Browse", command=browsefunc6).grid(row = 330, column = 3, sticky=W)
pathlabel = Label(top)

v2 = StringVar()
tkinter.Label(top, text = "= BLAST evalue 1").grid(row = 280, column = 3) 
tkinter.Entry(top, text=v2).grid(row = 280, column = 2) 
v2.set('"1e-10"')

v3 = StringVar()
tkinter.Label(top, text = "= BLAST evalue 3").grid(row = 320, column = 3) 
tkinter.Entry(top, text=v3).grid(row = 320, column = 2) 
v3.set('"1e-10"')

v28 = StringVar()
tkinter.Label(top, text = "= BLAST evalue 2").grid(row = 300, column = 3) 
tkinter.Entry(top, text=v28).grid(row = 300, column = 2) 
v28.set('"1e-10"')

tkinter.Label(top, text = "                ").grid(row = 0, column = 4)
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
tkinter.Label(top, text = "Expression Analysis:", fg="blue").grid(row = 0, column = 17, sticky=W)

tkinter.Checkbutton(top, text = "RSEM_blast1", variable = CV24, onvalue = 1, offvalue=0).grid(row = 20, column = 17, sticky=W)
tkinter.Label(top, text = "If Starting with Expression Analysis Provide Fasta Below:", fg='gray').grid(row = 10, column = 17, columnspan = 2) 
tkinter.Checkbutton(top, text = "Kallisto_blast1", variable = CV25, onvalue = 1, offvalue=0).grid(row = 30, column = 17, sticky=W)
tkinter.Checkbutton(top, text = "Salmon_blast1", variable = sal1, onvalue = 1, offvalue=0).grid(row = 31, column = 17, sticky=W)
v25 = StringVar()
tkinter.Entry(top, text=v25).grid(row = 20, column = 18) 
v25.set("")
browsebutton = Button(top, text="Browse", command=browsefunc8).grid(row = 20, column = 19, sticky=W)
pathlabel = Label(top)

v31 = IntVar()
tkinter.Checkbutton(top, text = "RSEM_blast2", variable = v31, onvalue = 1, offvalue=0).grid(row = 50, column = 17, sticky=W)
tkinter.Label(top, text = "If Starting with Expression Analysis Provide Fasta Below:", fg='gray').grid(row = 40, column = 17, columnspan = 2) 
v32 = IntVar()
tkinter.Checkbutton(top, text = "Kallisto_blast2", variable = v32, onvalue = 1, offvalue=0).grid(row = 60, column = 17, sticky=W)
tkinter.Checkbutton(top, text = "Salmon_blast2", variable = sal2, onvalue = 1, offvalue=0).grid(row = 61, column = 17, sticky=W)
v33 = StringVar()
tkinter.Entry(top, text=v33).grid(row = 50, column = 18) 
v33.set("")
browsebutton = Button(top, text="Browse", command=browsefunc11).grid(row = 50, column = 19, sticky=W)
pathlabel = Label(top)
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
tkinter.Label(top, text = "Remove Unnecessary Files:", fg="blue").grid(row = 300, column = 17, sticky=W)

tkinter.Checkbutton(top, text = "Cleanup", variable = CV28, onvalue = 1, offvalue=0).grid(row = 310,column = 17,sticky=W)
tkinter.Label(top, text = "Delete Intermediate Files").grid(row = 310, column = 18) 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
v16 = IntVar()
tkinter.Label(top, text = "Resources:", fg="blue").grid(row = 260, column = 17, sticky=W) 
tkinter.Label(top, text = "# of CPU Threads (0=Max) for Pincho").grid(row = 270, column = 17) 
tkinter.Entry(top, text=v16).grid(row = 270, column = 18) 
v16.set(0)

blast_resource = IntVar()
tkinter.Label(top, text = "# of CPU Threads (0=Max) for Blast").grid(row = 271, column = 17) 
tkinter.Entry(top, text=v16).grid(row = 271, column = 18) 
blast_resource.set(0)

v17 = IntVar()
tkinter.Label(top, text = "Max Memory [GB] (0=/=Max)", fg="red").grid(row = 280, column = 17) 
tkinter.Entry(top, text=v17).grid(row = 280, column = 18) 
v17.set(64)

tkinter.Label(top, text = "Make sure system swap space is compatible with max memory", fg="green").grid(row = 290, column = 17,  columnspan =2) 
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
button = tkinter.Button(text="                   START                   ", command=start).grid(row = 330, column = 17,sticky=W)
button = tkinter.Button(text="   High Quality Full Run   ", command=high_quality).grid(row = 330, column = 18,sticky=W)
button = tkinter.Button(text="  Test   ", command=test).grid(row = 320, column = 19,sticky=W)
button = tkinter.Button(text="  Clear  ", command=clear).grid(row = 330, column = 19,sticky=W)
button = tkinter.Button(text="QUIT", fg="red", command=quit).grid(row = 330, column = 20,sticky=W)
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
top.mainloop()


######################################################################################################################

######################################################################################################################
                                                ###   END August 16 2021
######################################################################################################################

######################################################################################################################