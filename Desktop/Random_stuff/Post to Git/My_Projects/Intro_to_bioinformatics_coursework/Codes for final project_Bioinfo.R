# Hb Project for Bioinformatics

# we have data for each Hb subunit (HBB1, HBB2, HBA1, HBAD), we have its sequence for 200 species of lizards, snakes, and turtles.
# we excluded data of crocodiles. each data for each subunit simply has the altitude it resides at, the consideration of whether 
# its high/low and the type of reptile it is.


# need results for the following:

# Amino acid consensus sequences of high living and low living of subunit HBB1 and HBB2, and its percent of amino acid residues

#RESULTS PART 1: will first do HBB1 of low living (to get consensus seq) and visualize the alignment 

## install DECIPHER
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DECIPHER")

##install ggmsa

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggmsa")

## make sure that DECIPHER and Biostrings are loaded
require(DECIPHER)
require(Biostrings)

#load your seq

HBB1_lowliving_Biostrings <- readDNAStringSet(filepath = "~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_low_living_snakes.fasta", 
                                       format = "fasta")

##alignment needs to be done 

HBB1_lowliving_aligned <- AlignSeqs(HBB1_lowliving_Biostrings)

##save alignment
writeXStringSet(HBB1_lowliving_aligned, file="HBB1_lowliving_aligned.fasta", format="fasta")

#visualize alignment

BrowseSeqs(HBB1_lowliving_aligned)

#save as PDF.

#visualize alignment via ggmsa

#HBB1_lowliving_aligned2 <- "~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_lowliving_aligned.fasta"

#my_lowliving_plot_1 <- ggmsa(msa = HBB1_lowliving_aligned2,
#                 start = 1,
#                  font = "helvetical",
#                color = "Chemistry_AA",
#           char_width = 0.9)

#my_lowliving_plot_1



## redo with high altitude living 

HBB1_highliving_Biostrings <- readDNAStringSet(filepath = "~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_high_living_snakes.fasta", 
                                       format = "fasta")

##alignment needs to be done 

HBB1_highliving_aligned <- AlignSeqs(HBB1_highliving_Biostrings)

##save alignment
writeXStringSet(HBB1_highliving_aligned, file="HBB1_highliving_aligned.fasta", format="fasta")

#visualize alignment

BrowseSeqs(HBB1_highliving_aligned)


##RESULTS PART 2: comparing consensus sequences:

#get the consensus sequences

consensus_lowliving <- consensusString(consensusMatrix(HBB1_lowliving_Biostrings))
consensus_highliving <- consensusString(consensusMatrix(HBB1_highliving_Biostrings))

#compare the sequences via pairwise alignment

pw_alignment <- pairwiseAlignment(consensus_lowliving, consensus_highliving, type="global")

#get alignment score
cat("Alignment Score:", pw_alignment@score, "\n")

#I got alignment score of 931.5754.

#alignment score is the number of matches between the aligned sequences 
#minus the number of mismatches and gaps. so 931 is fairly high. 


## translate Biostring DNAStringSet from nucleotide to amino acid 
#state vertebrate code for genetic code standard
vertebrate_code <- getGeneticCode("SGC0")

#translate high living
HBB1_highliving_aligned <- readDNAStringSet(filepath="~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_highliving_aligned.fasta", format="fasta")
HBB1_highliving_AA <- AlignTranslation(HBB1_highliving_aligned, type="AAStringSet", geneticCode = vertebrate_code)

BrowseSeqs(HBB1_highliving_AA)

writeXStringSet(HBB1_highliving_AA, file="HBB1_highliving_AA", format="fasta")

#translate low living 
HBB1_lowliving_aligned <- readDNAStringSet(filepath="~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_lowliving_aligned.fasta", format="fasta")
HBB1_lowliving_AA <- AlignTranslation(HBB1_lowliving_aligned, type="AAStringSet", geneticCode = vertebrate_code)

BrowseSeqs(HBB1_lowliving_AA)

writeXStringSet(HBB1_lowliving_AA, file="HBB1_lowliving_AA.fasta", format="fasta")


##compare only crotalus high vs low altitude 
HBB1_Crotalus_low_AA <- readAAStringSet(filepath="~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_Crotalus_low_AA.fasta", format="fasta")
BrowseSeqs(HBB1_Crotalus_low_AA)

HBB1_Crotalus_high_AA <- readAAStringSet(filepath="~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_Crotalus_high_AA.fasta", format="fasta")
BrowseSeqs(HBB1_Crotalus_high_AA)


###compare GC content

library(seqinr)

#HBB1_lowliving_Biostrings <- readDNAStringSet(filepath = "~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_low_living_snakes.fasta", 
                                       #format = "fasta")

input_HBB1_low <- "~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_low_living_snakes.fasta"

GCinput<- read.fasta(input_HBB1_low,as.string =TRUE, seqtype= "DNA", seqonly=FALSE)
#empty dataframe to hold results
dataframeGCout <- data.frame()

#work through each sequence 
for(k in 1: length(GCinput)){
	sequencetoGC <-getSequence(GCinput[k])
	sequencetoGC <- unlist (sequencetoGC)
	GCout <- GC(sequencetoGC)
	dataframeGCout[k,1]<- getName (GCinput[k])
	dataframeGCout[k,2]<-GCout}


#Print results
colnames(dataframeGCOut)<- c("Sequence", "GC Content")

dataframeGCout

dataframeGCout <- write.xlsx(dataframeGCout,
                                         sheetName = "GCout_sheet",
                                              file = "dataframeGCout.xlsx")

#HBB1_high
input_HBB1_high <- "~/Desktop/Research Projects/Hb-O2_altitude_project/HB_Project/HBB1_high_living_snakes.fasta"

GCinput<- read.fasta(input_HBB1_high,as.string =TRUE, seqtype= "DNA", seqonly=FALSE)
#empty dataframe to hold results
dataframeGCout <- data.frame()

#work through each sequence 
for(k in 1: length(GCinput)){
	sequencetoGC <-getSequence(GCinput[k])
	sequencetoGC <- unlist (sequencetoGC)
	GCout <- GC(sequencetoGC)
	dataframeGCout[k,1]<- getName (GCinput[k])
	dataframeGCout[k,2]<-GCout}


#Print results
colnames(dataframeGCout)<- c("Sequence", "GC Content")

dataframeGCout

#put into excel
dataframeGCout <- write.xlsx(dataframeGCout,
                                         sheetName = "GCout_sheet",
                                              file = "dataframeGCout.xlsx")
