#!/usr/bin/env Rscript
library("dplyr")
#library("data.table")
library("ACE")
library("QDNAseq.hg38")
library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)
patient <- args[1]

# find ace rds output
my_rds <- list.files(".", recursive = TRUE, pattern = "*1000kbp.rds$")

# read rds and name after sample
object_name <- sub("_filter_90_150/1000kbp.rds", "", my_rds)
for (i in seq_along(my_rds)) {
  assign(paste0(patient,"_",object_name[i]), readRDS(my_rds[i]))
}

# run getadjustedsegments and get output, name as $sample
for (my_object in object_name) {
  current_object <- get(paste0(patient,"_",my_object))
  model1 <- singlemodel(current_object, QDNAseqobjectsample = 1)
  bestfit1 <- model1$minima[tail(which(model1$rerror == min(model1$rerror)), 1)]
  assign(paste0(patient,"_", my_object, "_segments"), 
                getadjustedsegments(current_object,
                                   QDNAseqobjectsample = 1,
                                   cellularity = bestfit1 ))
}

# add sample_name column and name as $sample_segments
for (my_object in object_name) {
  current_object <- get(paste0(patient, "_", my_object, "_segments"))
  current_object$sample_id <- paste0(patient, "_", my_object)
  assign(paste0(patient, "_", my_object, "_segments"), current_object, envir = .GlobalEnv)
}

# this is necessary to expand all the segments to 1Mbp and then output to csv
options(scipen = 999)
step_size <- 1000000
final_df <- data.frame("sample_id" = c(), "Chromosome" = c(), "Start" = c(), "End" = c(), "Copies" = c())
for (sample in (ls(pattern = "segments"))) {
  SampleID <- sub("_segments$", "", sample)
  sample <- get(sample)
  Chromosomes <- unique(sample$Chromosome)
  chromosome_df <- data.frame("sample_id" = c(), "Chromosome" = c(), "Start" = c(), "End" = c(), "Copies" = c())
  for ( chr in Chromosomes ) {
    df <- sample[sample$Chromosome == chr,] %>%
      select(Chromosome,Start,End,Copies) %>%
      mutate(Combined = paste(Start, End, Copies, sep = "_")) %>%
      group_by(Combined) %>%
      mutate(group_indices = cur_group_id()) %>%
      ungroup()
   
    split_df_list <- split(df, df$group_indices)

    result_df <- data.frame("sample_id" = c(), "Chromosome" = c(), "Start" = c(), "End" = c(), "Copies" = c())
    for (i in seq_along(split_df_list)) {
      num_intervals <- ceiling((split_df_list[[i]]$End - split_df_list[[i]]$Start + 1 + 1e-6) / step_size)
      seq_values <- seq(split_df_list[[i]]$Start, by = step_size, length.out = num_intervals)
      stops <-ifelse(seq_values + 999999 < max(split_df_list[[i]]$End),seq_values+999999,max(split_df_list[[i]]$End))
      processed_df <- split_df_list[[i]] %>% uncount(weights = num_intervals) %>%
      mutate(sample_id = SampleID, Start = seq_values, End= stops) 
      result_df <- bind_rows(result_df, processed_df)
    }
  chromosome_df <- bind_rows(chromosome_df, result_df)
  }
  final_df <- bind_rows(final_df, chromosome_df)
}

# keep only the chrom start end combinations that are in all the samples
filtered_df <- final_df %>% group_by(Chromosome,Start,End) %>% filter(n() >= length(my_rds))

# rename and select columns for medicc2 
out_df <- filtered_df %>% 
  rename(chrom = Chromosome, start = Start, end = End) %>% 
  select(sample_id, chrom, start, end, Copies)

write.table(out_df,file=paste0(patient,".tsv"),quote = FALSE,
            sep = "\t", row.names = FALSE)
