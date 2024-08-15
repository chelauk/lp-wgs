#!/usr/bin/env Rscript
library("dplyr")
#library("data.table")
library("ACE")
library("QDNAseq.hg38")
library("tidyverse")
options(scipen = 999) # prevent scientific notation

## functions

explode_ranges <- function(df, step_size) {
  expanded_df <- data.frame()
  
  for (i in 1:nrow(df)) {
    
    Chromosome <- df$Chromosome[i]
    Start <- df$Start[i]
    End <- df$End[i]
    Copies <- df$Copies[i]
    
    seq_Starts <- seq(Start, End - step_size, by = step_size)
    seq_Ends <- seq(Start + step_size - 1, End, by = step_size)
    if ( seq_Ends[length(seq_Ends)] < End ) {
      seq_Ends <- c(seq_Ends, End)
      print(paste("last two seq_Ends", seq_Ends[length(seq_Ends) - 1], seq_Ends[length(seq_Ends)]))
    }
    if ( seq_Starts[length(seq_Starts)] + step_size < End){
      seq_Starts <- c(seq_Starts,seq_Starts[length(seq_Starts)] + step_size)
    }
    
    Chromosome <- rep(Chromosome,length(seq_Starts))
    Copies <- rep(Copies,length(seq_Starts))
    #print(paste("Chromosome: ", Chromosome[1]))
    #print(paste("length chromosome: ", length(Chromosome)))
    #print(paste("first start:",seq_Starts[1], "last start:", seq_Starts[length(seq_Starts)]))
    #print(paste("length seq_Starts: ",length(seq_Starts)))
    #print(seq_Starts)
    #print(seq_Ends)
    #print(paste("first end:",seq_Ends[1], "last end:", seq_Ends[length(seq_Ends)]))
    #print(paste("length seq_Ends: ",length(seq_Ends)))
    #print(paste("length Copies: ",length(Copies)))
    new_rows <- data.frame(
      Chromosome = Chromosome,
      Start = seq_Starts,
      End = seq_Ends,
      Copies = Copies
    )
    
    expanded_df <- bind_rows(expanded_df, new_rows)
  }
  
  return(expanded_df)
}

###

args <- commandArgs(trailingOnly = TRUE)
patient <- args[1]

# find ace rds output
my_rds <- list.files(".", recursive = TRUE, pattern = "*1000kbp.rds$")

# this will give us a vector of sample_names
object_name <- sub("_filter_.*\\/1000kbp\\.rds", "", my_rds)
sample_list <- c()
for (i in seq_along(my_rds)) {
  sample_list <- c(sample_list, paste0(patient,"_", object_name[i]))
  assign(paste0(patient,"_",object_name[i]), readRDS(my_rds[i]))
}

# run get adjustedsegments and get output, name as $sample
# medicc2 parameter info will be stored
medicc2_params <- "medicc2_params.txt"
for (my_object in object_name) {
  # get ace output and best fits using minima
  current_object <- get(paste0(patient,"_",my_object))
  model1 <- singlemodel(current_object, QDNAseqobjectsample = 1, ploidy = 2 )
  # model1 is an object with different minima, trouble is there are 
  # options either bestfit or lastminima we need to be able to adjust ploidy
  bestfit1 <- model1$minima[tail(which(model1$rerror == min(model1$rerror)), 1)]
  # I need to work on this
  log_txt <- paste0(patient, "_", my_object, " cellularity: ", bestfit1, ", ploidy: ", model1$ploidy, "\n")
  cat(log_txt, file = medicc2_params, append = TRUE)
  assign(paste0(patient,"_", my_object, "_segments"), 
                getadjustedsegments(current_object,
                                   QDNAseqobjectsample = 1,
                                   cellularity = bestfit1 ))
}



# add sample_name column and name as $sample_segments
for (my_object in object_name) {
  current_object <- get(paste0(patient, "_", my_object, "_segments"))
  #current_object <- segmentstotemplate(current_object, meanci = 8)[,c(2,3,4,5)]
  current_object <- explode_ranges(current_object, 1000000)
  all_bins <- get(paste0(patient, "_", my_object))
  all_bins  <- as.data.frame(all_bins@assayData$copynumber)
  all_bins <- all_bins %>% rownames_to_column(var = "position") %>%
    select(position) %>%
    separate(position, into = c("Chromosome", "Start","End"), sep = "[:-]") %>%
    mutate( Start = as.numeric(Start), End = as.numeric(End))
  
  current_object <- left_join(all_bins,current_object,by=c("Chromosome","Start","End"))
  #current_object <- mutate(Chromosome = as.numeric(Chromosome))
  current_object$sample_id <- paste0(patient, "_", my_object)
  print(head(current_object))
  current_object <- current_object %>% dplyr::select("sample_id","Chromosome","Start","End","Copies") %>% 
    rename(chrom = Chromosome, start = Start, end = End)
  assign(paste0(patient, "_", my_object, "_segments"), current_object, envir = .GlobalEnv)
}
dataframes <- mget(ls(pattern = "segments"))
out_df <- bind_rows(dataframes)
out_df <- pivot_wider(out_df,names_from = "sample_id", values_from = "Copies")
out_df <- out_df[complete.cases(out_df),]
my_names<- paste0(patient,"_",object_name)

out_df <- out_df %>% group_by(chrom, across(my_names)) %>%
  summarize(
    start = min(start),
    end = max(end)
  )
out_df <- out_df %>%
  pivot_longer(
    cols = -c(chrom,start,end),    # Specify columns to gather
    names_to = "sample_id",      # Name for the new key column
    values_to = "Copies"    # Name for the new value column
  )
out_df$Diploid <- 2
out_df <- mutate(out_df, chrom = as.integer(chrom))
out_df <- out_df %>% arrange(chrom,start)
#raise MEDICCIOError("The samples have different segments!\n"
#medicc.io.MEDICCIOError: The samples have different segments!
#Total number of unique segments: 3356
write.table(out_df,file=paste0(patient,".tsv"),quote = FALSE,
            sep = "\t", row.names = FALSE)
