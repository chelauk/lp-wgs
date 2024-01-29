#!/usr/bin/env Rscript
library("dplyr")
#library("data.table")
library("ACE")
library("QDNAseq.hg38")

args <- commandArgs(trailingOnly = TRUE)
patient <- args[1]

# find ace rds output
my_rds <- list.files(".", recursive = TRUE, pattern = "*1000kbp.rds$")

# read rds and name after sample
object_name <- sub("/low_pass_wgs/ace_1000kb/filter_90_150/1000kbp.rds",
                   "", my_rds)
for (i in seq_along(my_rds)) {
  assign(object_name[i], readRDS(my_rds[i]))
}

# run getadjustedsegments and get output, name as $sample
for (my_object in object_name) {
  current_object <- get(my_object)
  model1 <- singlemodel(current_object, QDNAseqobjectsample = 1)
  bestfit1 <- model1$minima[tail(which(model1$rerror == min(model1$rerror)), 1)]
  assign(paste0(my_object, "_segments"), getadjustedsegments(current_object,
                                        QDNAseqobjectsample = 1,
                                        cellularity = bestfit1 ))
}

# add sample_name column and name as $sample_segments
for (my_object in object_name) {
  current_object <- get(paste0(my_object, "_segments"))
  current_object$sample_id <- my_object
  assign(paste0(my_object, "_segments"), current_object, envir = .GlobalEnv)
}

# concatenate all segment dfs
my_dfs <- my_dfs <- mget(ls(pattern = "segments"))
combined_df <- dplyr::bind_rows(my_dfs)

# Group by Chromosome, Start, and End and count occurrences
occurrence_counts <- combined_df %>%
  group_by(Chromosome, Start, End) %>%
  summarise(Count = n())

# keep only the chrom start end combinations that are in all the dfs
filtered_df <- combined_df %>%
  inner_join(occurrence_counts, by = c("Chromosome", "Start", "End")) %>%
  filter(Count == length(my_rds))

# rename and select columns for medicc2 
out_df <- filtered_df %>% 
  rename(chrom = Chromosome, start = Start, end = End) %>% 
  select(sample_id, chrom, start, end, Copies)

write.table(out_df,file=paste0(patient,".tsv"),quote = FALSE,
            sep = "\t", row.names = FALSE)
