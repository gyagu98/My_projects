#reconstructing mtgenome of planaria

#make folder for the reads, and inside, put
	#the binaries of minimap2 and samtools
	#raw reads
	#reference seed of species that is most closely related (found from BLAST in NCBI)

cd minimap2

#minimap codes, make 1 for the forward read and 1 for the reverse read (raw files). these codes will make the map 
#against the reference

#here, is the reference against the entire mitochondrial genome
##edit the codes in according to the F number, also, run these separately 

./minimap2 -a ../amaga_mito_ref.fasta ../F525_R1_001.fastq.gz   -axsr > ../F525_R1_001.sam

./minimap2 -a ../amaga_mito_ref.fasta ../F525_R2_001.fastq.gz  -axsr > ../F525_R2_001.sam

#here is the reference against the cytochrome gene

./minimap2 -a ../amaga_cyto_ref.fasta ../F525_R1_001.fastq.gz -axsr > ../F525_R1_001.sam

./minimap2 -a ../amaga_cyto_ref.fasta ../F525_R2_001.fastq.gz -axsr > ../F525_R2_001.sam

#go back directory
cd ..

#now, we will use samtools to clean/filter likely for size or length

#these two codes first convert between sam and bam files
samtools view -b -F 4 F525_R1_001.sam >  F525_R1_mapped.bam

samtools view -b -F 4 F525_R2_001.sam >  F525_R2_mapped.bam

#then, these next two codes will convert to fastq files
samtools fastq F525_R1_mapped.bam > F525_R1_filter.fastq

samtools fastq F525_R2_mapped.bam > F525_R2_filter.fastq

#those filter.fastq files are needed because we will then use them to run Pincho
#create a folder with the F number and the R1 and R2 reads separate, and in that folder, open Pincho

'/home/tuvok/Documents/Pincho_v01_Master/Pincho_v02_conf.py' 

#we will run Pincho the following way:
	#1 hit high quality run
	#2 remove trimmomatic, Rcorrector, Rnaspades
	#3 click choose file of only forward run (we do this bc R1 and R2 are different lengths; for R2 also choose forward read)
	#4 click run. the first run will compress the file will take ~30 seconds. Remove the folder made.
	#5 then for the second run, you will choose the new compressed file, it should take  ~5 min to complete, and when it is done, you want the entire folder of pincho results