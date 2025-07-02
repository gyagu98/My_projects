#planaria alignments

require(DECIPHER)
require(Biostrings)

setwd("~/Desktop/Research Projects/planaria")

Doroceras_Biostrings_set <- readDNAStringSet(filepath = "~/Desktop/Research Projects/planaria/Deroceras_phylogeny/Deroceras_COX1_collection.fasta", 
                                       format = "fasta")

Doroceras_aligned <- AlignSeqs(Doroceras_Biostrings_set)


BrowseSeqs(Doroceras_aligned)

writeXStringSet(Doroceras_aligned, file="Dorocerase_aligned.fasta", format="fasta")

##repeat for Amaga

Amaga_Biostrings_set <- readDNAStringSet(filepath = "~/Desktop/Research Projects/planaria/Amaga_phylogeny/amaga_COX_UBRP_V2.fasta", 
                                       format = "fasta")

Amaga_aligned <- AlignSeqs(Amaga_Biostrings_set)


BrowseSeqs(Amaga_aligned)

writeXStringSet(Amaga_aligned, file="Amaga_aligned.fasta", format="fasta")


#########################run iqtree2:

require(doParallel)
require(ape)

 iqtree_runner_R <- function(input_alignment_file_user,
                           iqtree_bin_location_user,
                           run_iqtree_call_user = c("serial","parallel"),
                           iqtree_quiet_user = FALSE,
                           iqtree_redo_user = FALSE,
                           iqtree_prefix_user = NULL,
                           iqtree_input_model_m_user = NULL,
                           iqtree_mtree_user = FALSE,
                           iqtree_sequence_type_st_user = NULL,
                           iqtree_UFBoot_B_user = NULL,
                           iqtree_UFBoot_bnni_user = TRUE,
                           iqtree_nonpara_boot_b_user = NULL,
                           iqtree_SHaLRT_alrt_user = NULL,
                           iqtree_multicore_T_user = NULL,
                           iqtree_maxcores_ntmax_user = NULL,
                           iqtree_partition_user = NULL,
                           iqtree_partition_if_by_position_user = TRUE,
                           iqtree_partition_file_user = NULL,
                           iqtree_chronogram_from_tree_file_user = NULL,
                           iqtree_chronogram_mrca_bounds_user = NULL,
                           iqtree_date_file_user = NULL,
                           iqtree_pass_parameters_user = NULL){

## require libraries

require(doParallel)
require(ape)

############# input from user

## sequence and tree arguments

input_files <- input_alignment_file_user 

# iqtree location

iqtree_bin_location <- iqtree_bin_location_user

# serial or parallel calls

run_iqtree <- run_iqtree_call_user

iqtree_quiet <- iqtree_quiet_user

####### basic iqtree parameters

# overwrite all previous output: -redo

iqtree_redo <- iqtree_redo_user

# add prefi: --prefix

iqtree_prefix <- iqtree_prefix_user

# add model: -m

iqtree_model <- iqtree_input_model_m_user

# increase accuracy: -mtree

iqtree_mtree <- iqtree_mtree_user

# Specify sequence type as either of NULL (autodetect) DNA, AA, BIN, MORPH, CODON or NT2AA for DNA

iqtree_st <- iqtree_sequence_type_st_user

# UFBoot mode: -B 1000

iqtree_UFBoot_B <- iqtree_UFBoot_B_user

# Red_impact_UFBoot

iqtree_UFBoot_bnni <- iqtree_UFBoot_bnni_user

# nonparametric bootstrap -b 100

iqtree_nonpara_b <- iqtree_nonpara_boot_b_user

# Assessing branch supports with single branch tests SH-like approximate likelihood ratio test: -alrt 1000

iqtree_alrt <- iqtree_SHaLRT_alrt_user

# multicore -T AUTO -T 2

iqtree_multicore_T <- iqtree_multicore_T_user

# maxcores -ntmax 

iqtree_maxcores_ntmax <- iqtree_maxcores_ntmax_user

##############       iqtree_partition_user

iqtree_partition <- iqtree_partition_user
do_by_codon_pos <- iqtree_partition_if_by_position_user


if(!is.null(iqtree_partition)) {

      iqtree_partition_matrix <- do.call(rbind,iqtree_partition)
             partitions_names <- rownames(iqtree_partition_matrix)
              type_partitions <- iqtree_partition_matrix[,1]
                   start_part <- as.numeric(iqtree_partition_matrix[,2])
                     end_part <- as.numeric(iqtree_partition_matrix[,3])
                   codon_part <- iqtree_partition_matrix[,4]
          iqtree_partition_df <- data.frame(partition_name = partitions_names, type = type_partitions, start = start_part, end = end_part, codon = codon_part, stringsAsFactors = FALSE)
rownames(iqtree_partition_df) <- NULL

# output file name

out_part_file <- sub(".*/", "", input_files)
out_part_file <- sub("[.]", "_", out_part_file)
out_part_file <- paste0(out_part_file,"_partition.txt")


for(i in 1:nrow(iqtree_partition_df)) {
  # i <- 3
   one_line <- paste0(iqtree_partition_df$type[i], ", ", iqtree_partition_df$partition_name[i], " = ", 
                      iqtree_partition_df$start[i], "-", iqtree_partition_df$end[i])

if(!is.na(iqtree_partition_df$codon[i])) {

# do_by_codon_pos
   if(do_by_codon_pos) {
              pos1_start <- iqtree_partition_df$start[i]
               pos1_name <- paste0(iqtree_partition_df$partition_name[i],"_pos1")
              pos2_start <- iqtree_partition_df$start[i]+1
               pos2_name <- paste0(iqtree_partition_df$partition_name[i],"_pos2")
              pos3_start <- iqtree_partition_df$start[i]+2
               pos3_name <- paste0(iqtree_partition_df$partition_name[i],"_pos3")

   one_line <- paste0(iqtree_partition_df$type[i], ", ", pos1_name, " = ", pos1_start, "-", iqtree_partition_df$end[i], "\\", iqtree_partition_df$codon[i], "\n",
                      iqtree_partition_df$type[i], ", ", pos2_name, " = ", pos2_start, "-", iqtree_partition_df$end[i], "\\", iqtree_partition_df$codon[i], "\n",
                      iqtree_partition_df$type[i], ", ", pos3_name, " = ", pos3_start, "-", iqtree_partition_df$end[i], "\\", iqtree_partition_df$codon[i], "\n")
                       }

                                        }

#     one_line <- paste0(one_line, "\\", iqtree_partition_df$codon[i]) }

cat(one_line,file=out_part_file,sep="\n",append=TRUE)

                                      }
                                      rm(i)

iqtree_partition_files <- out_part_file 

                                         } else {

iqtree_partition_files <- iqtree_partition_file_user

                                         }

##############       chronogram_functions

iqtree_chronogram_from_tree <- iqtree_chronogram_from_tree_file_user

iqtree_mrca_bounds <- iqtree_chronogram_mrca_bounds_user

if(!is.null(iqtree_mrca_bounds)) {

# name for DATE_FILE

out_part_file <- sub(".*/", "", input_files)
out_part_file <- sub("[.]", "_", out_part_file)
out_date_file <- paste0(out_part_file,"_dates.txt")

# build DATE_FILE

n_total <- length(iqtree_mrca_bounds)
n_constraints <- seq(1,by=2, to=n_total)

for(i in n_constraints) {

# i <- 1

   taxa_part <- paste0(iqtree_mrca_bounds[[i]], collapse=",")
   date_part <- paste0("-", iqtree_mrca_bounds[[i+1]])

   one_line <- paste0(taxa_part, " ", date_part)

#     one_line <- paste0(one_line, "\\", iqtree_partition_df$codon[i]) }

cat(one_line,file=out_date_file,sep="\n",append=TRUE)

                                      }
                                      rm(i)

iqtree_date_files <- out_date_file 

                                 } else {

iqtree_date_files <- iqtree_date_file_user

                                         }

# iq pass paramters

iqtree_pass <- iqtree_pass_parameters_user

# get work directory to go back

master_wd_directory <- getwd()


## get iqtree path

path_of_iqtree <- paste0(iqtree_bin_location)

######     prepare iqtree commands

list_iqtree_scripts <- list()
name_of_analyses_iqtree <- character()

counter <- 0


   for(i in 1:length(input_files)) {
      # i <- 1
    one_alignment_file <- input_files[i]
    counter <- counter + 1

# simplify names

one_input_file <- sub(".*/", "", one_alignment_file)

# prepare commands

line_iqtree <- paste0("iqtree2  -s ", one_input_file, "  ")

# list of iqtre parameters

list_of_iqtree <- c(line_iqtree)

counter_2 <- 1

# Specify sequence type as either of NULL (autodetect) DNA, AA, BIN, MORPH, CODON or NT2AA for DNA

if(!is.null(iqtree_st)) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  -st  ", iqtree_st, "  ") 
   list_of_iqtree[counter_2] <- out_line
                     }

# add prefi: --prefix

if(!is.null(iqtree_prefix)) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  --prefix  ", iqtree_prefix, "  ") 
   list_of_iqtree[counter_2] <- out_line
                     }

# add model: -m

if(!is.null(iqtree_model)) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  -m  ", iqtree_model, "  ") 
   list_of_iqtree[counter_2] <- out_line
                     }

# increase accuracy: -mtree

if(iqtree_mtree) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  -mtree  ") 
   list_of_iqtree[counter_2] <- out_line
                     }

# overwrite all previous output: -redo

if(iqtree_redo) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  -redo  ") 
   list_of_iqtree[counter_2] <- out_line
                     }

# add partition

if(!is.null(iqtree_partition_files)) {counter_2 <- counter_2 + 1
                   out_line <- paste0("  -p  ", iqtree_partition_files[i]) 
   list_of_iqtree[counter_2] <- out_line

                              }

# add DATE_FILE

if(is.null(iqtree_chronogram_from_tree)) {

if(!is.null(iqtree_date_files)) {counter_2 <- counter_2 + 1
                   out_line <- paste0("  --date  ", iqtree_date_files[i], ' --date-tip 0') 
   list_of_iqtree[counter_2] <- out_line

                              }

                                        }

# make chronogram from tree

if(!is.null(iqtree_chronogram_from_tree)) {counter_2 <- counter_2 + 1
                   out_line <- paste0("  --date  ", iqtree_date_files[i], ' --date-tip 0', " -te ", iqtree_chronogram_from_tree) 
   list_of_iqtree[counter_2] <- out_line

                              }

##### ultra fast bootstrap

# UFBoot mode: -B 1000

if(!is.null(iqtree_UFBoot_B)) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  -B ", iqtree_UFBoot_B, "  ") 
   list_of_iqtree[counter_2] <- out_line
                     }

# Red_impact_UFBoot

if(!is.null(iqtree_UFBoot_B)) {
         if(iqtree_UFBoot_bnni) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  -bnni  ") 
   list_of_iqtree[counter_2] <- out_line
                                }
                              }

##### non-parameteric bootstrap

# nonparametric bootstrap -b 100

if(!is.null(iqtree_nonpara_b)) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  -b ", iqtree_nonpara_b, "  ") 
   list_of_iqtree[counter_2] <- out_line
                     }

##### SH-like approximate likelihood ratio test

# Assessing branch supports with single branch tests : -alrt 1000

if(!is.null(iqtree_alrt)) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  -alrt ", iqtree_alrt, "  ") 
   list_of_iqtree[counter_2] <- out_line
                     }

### multicore

# multicore -T AUTO -T 2

if(!is.null(iqtree_multicore_T)) { counter_2 <- counter_2 + 1
                   out_line <- paste0("  -T ", iqtree_multicore_T, "  ") 
   list_of_iqtree[counter_2] <- out_line
                     }

if(!is.null(iqtree_multicore_T)) {
          if(!is.null(iqtree_maxcores_ntmax)) {
                                    counter_2 <- counter_2 + 1
                                     out_line <- paste0("  -ntmax ", iqtree_maxcores_ntmax, "  ") 
                    list_of_iqtree[counter_2] <- out_line
                                              }
                                  }


### iqtree pass paramters

if(!is.null(iqtree_pass)) { counter_2 <- counter_2 + 1
                   out_line <- paste0(" ", iqtree_pass, "  ") 
   list_of_iqtree[counter_2] <- out_line
                          }


# quiet
                        
if(iqtree_quiet) { counter_2 <- counter_2 + 1
                   out_line <- paste0(" -quiet  ") 
   list_of_iqtree[counter_2] <- out_line
                          }

### all lines together

iqtree_out <- paste0(list_of_iqtree, collapse = " ")
print(iqtree_out)

  list_iqtree_scripts[[counter]] <- iqtree_out
name_of_analyses_iqtree[counter] <- one_input_file

                                  }
                         rm(i)
                       if(exists("counter")){rm(counter)}
                       if(exists("counter_2")){rm(counter_2)}


#####################################################################################
#### open and run iqtree from R

# create vector to export path -- export PATH=/opt/local/bin:/opt/local/sbin:$PATH

export_path_iqtree_bins <- normalizePath(iqtree_bin_location)
export_path_iqtree_vector <- paste0("export PATH=",export_path_iqtree_bins,":$PATH")

## iqtree execute function

execute_iqtree_terminal <- function(export_path_section = export_path_iqtree_vector, 
                                   iqtree_command, 
                                   intern = FALSE, 
                                    wait = FALSE){

system(paste0(export_path_section, "\n",
                    iqtree_command, "\n"), 
                                intern = intern, 
                                  wait = wait)
                                                 }

## run list elements with iqtree commands

## serial
# run_iqtree <- "serial"

if(run_iqtree == "serial") {

cat("\n****** started serial analisis of: ", length(list_iqtree_scripts), "\n" )

for(j in 1:length(list_iqtree_scripts)) {

cat("***** iqtree analysis name: ",name_of_analyses_iqtree[j], "\n")
cat("***** iqtree commands: ", list_iqtree_scripts[[j]], "\n")

execute_iqtree_terminal(export_path_section = export_path_iqtree_vector, 
                             iqtree_command =  list_iqtree_scripts[[j]], 
                                    intern = FALSE, 
                                      wait = TRUE)
cat("***** DONE ***** \n")
                                        }
                                        rm(j)

                           }


### parallel
# run_iqtree <- "parallel"

if(run_iqtree == "parallel") {

cat("\n****** started parallel analisis of: ", length(list_iqtree_scripts), "\n" )
print(name_of_analyses_iqtree)

ncores <- parallel::detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)
foreach(iqtree_i=list_iqtree_scripts, .export = "execute_iqtree_terminal") %dopar% 
                               {execute_iqtree_terminal(iqtree_command=iqtree_i)}
parallel::stopCluster(cl)

cat("\n ****** DONE ****** \n")

                           }

### to return to user grab treefile

iqtrees_vector <- list.files(pattern = "\\.treefile$", ignore.case=TRUE)

list_of_trees <- list()

for(i in 1:length(iqtrees_vector)) {list_of_trees[[i]] <- ape::read.tree(iqtrees_vector[i])}
rm(i)

### if chronograms

iqtree_chronograms <- list.files(pattern = "\\.timetree.nex$", ignore.case=TRUE)

if(length(iqtree_chronograms) > 0) {

for(i in 1:length(iqtree_chronograms)) {
                               one_chronogram <- ape::read.nexus(iqtree_chronograms[i])
                                      n_trees <- length(list_of_trees)
                   list_of_trees[[n_trees+1]] <- one_chronogram }
rm(i)

                                   }

############ open and run iqtree from R: DONE

setwd(master_wd_directory)

cat("\n\n *************         iqtree outfiles written to        ************* \n")
print(master_wd_directory)

return(list_of_trees)

                                           }
##    END OF FUNCTION   


####### TO RUN: SIMPLE iqtree2 run (deroceras phylogeny)

setwd("~/Desktop/Research Projects/planaria/Deroceras_phylogeny")

my_path_to_iqtree2 <- "~/Desktop/my_programs"


COI_simple_phylogeny <- iqtree_runner_R (input_alignment_file_user = "~/Desktop/Research Projects/planaria/Deroceras_phylogeny/Dorocerase_aligned.fasta",
                           iqtree_bin_location_user = my_path_to_iqtree2,
                           run_iqtree_call_user = "serial",
                           iqtree_input_model_m_user = "MFP+MERGE")


## iqtree2 run with defined codon positions and get 1000 bootstraps (see the best fit model in the console "corrected akaike ...")
COI_with_partitions_phylogeny_boot <- iqtree_runner_R (input_alignment_file_user = "~/Desktop/Research Projects/planaria/Deroceras_phylogeny/Dorocerase_aligned.fasta",
                           iqtree_bin_location_user = my_path_to_iqtree2,
                           run_iqtree_call_user = "serial",
                           iqtree_input_model_m_user = "TVM+F+I+I+R2",
                           iqtree_UFBoot_B_user = 1000,
                           iqtree_UFBoot_bnni_user = TRUE)

###can visualize .treefile via itol.
#reroot the tree on limax, light grey dashes are, read about what is limax? (simple tree)

##visualize boot strap .ufboot via


######redo with amaga:

setwd("~/Desktop/Research Projects/planaria/Amaga_phylogeny")

my_path_to_iqtree2 <- "~/Desktop/my_programs"


COI_simple_phylogeny <- iqtree_runner_R (input_alignment_file_user = "~/Desktop/Research Projects/planaria/Deroceras_phylogeny/Amaga_aligned.fasta",
                           iqtree_bin_location_user = my_path_to_iqtree2,
                           run_iqtree_call_user = "serial",
                           iqtree_input_model_m_user = "MFP+MERGE")


## iqtree2 run with defined codon positions and get 1000 bootstraps (see the best fit model in the console "corrected akaike ...")
COI_with_partitions_phylogeny_boot <- iqtree_runner_R (input_alignment_file_user = "~/Desktop/Research Projects/planaria/Deroceras_phylogeny/Amaga_aligned.fasta",
                           iqtree_bin_location_user = my_path_to_iqtree2,
                           run_iqtree_call_user = "serial",
                           iqtree_input_model_m_user = "GTR+F+I+G4",
                           iqtree_UFBoot_B_user = 1000,
                           iqtree_UFBoot_bnni_user = TRUE)


















