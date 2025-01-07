##TO RUN AND VISUALIZE THE MT GENOME RECONSTRUCTION:: 

#We will be doing mitoz. In the folder, we need 3 files: R1 and R2 raw reads, and the fasta file of the reconstruction mtDNA
#but, before we run we have to make sure the fasta file has to:
	# - have no insertions
	# - short fasta name (10 characters minimum)
	# - keep 1st sequence of the fasta file only, and remove the others. make the file into 1 sequence

#open terminal in the folder that has those three files 

conda activate mitozEnv 

# the part that says clade, we choose the clade most close to the animal. planaria's phylum is Platyhelminthes, and the closest to that is either Annelida (segmented worms (this appears closer bc it has the same clade in the phylogeny as Phatyhelmintes)) 
#or Nematoda (roundworms). for the slug, whose phylum is Mollusca, we will choose that bc it is already an option in the code of mitoz 

#run each code within each folder: 


mitoz annotate --outprefix Gretta_mt --thread_number 12 --fastafiles Amaga_sp_F523_mtDNA_reconstruction.fasta --fq1 F523_R1_001.fastq.gz --fq2 F523_R2_001.fastq.gz --species_name planaria_1 --genetic_code 9 --clade Annelida-segmented-worms

mitoz annotate --outprefix Gretta_mt --thread_number 12 --fastafiles Amaga_sp_F525_mtDNA_reconstruction.fasta --fq1 F525_R1_001.fastq.gz --fq2 F525_R2_001.fastq.gz --species_name planaria_1 --genetic_code 9 --clade Annelida-segmented-worms

#the slug has no raw reads because we did not get it sequenced
mitoz annotate --outprefix Gretta_mt --thread_number 12 --fastafiles deroceras_reticulatum_mtDNA.fasta --species_name doceras --genetic_code 9 --clade Mollusca

##understanding the output:
# the summary.txt tells you the genes found, their locations, and the potentially missing genes. CDS (coding region), tRNA (transfer RNA region), rRNA (ribosomal RNA region)
# we can also go into the mt_visualization directory and get the png of the mt reconstruction. around the genome are the genes, and there is a frequency of how much the gene was found. 
# the fasta reconstruction mitoscaf.fa.gbf file is the one that is used to submit to NCBI upon certain edits