#!/usr/bin/env Rscript
library("dplyr")
library("ACE")
library("QDNAseq.hg38")
library("tidyverse")
options(scipen = 999) # prevent scientific notation

args <- commandArgs(trailingOnly = TRUE)
patient <- args[1]
input_ploidys <- as.integer(unlist(strsplit(args[2]," ")))
bin_dir <- args[3]

source(paste0(bin_dir,"/00_general_functions.R")) # import functions
# find ace rds output
my_rds <- list.files(".", recursive = TRUE, pattern = "*1000kbp.rds$")

# this will give us a vector of sample_names
object_name <- sub("_filter_.*\\/1000kbp\\.rds", "", my_rds)
sample_list <- c()
for (i in seq_along(my_rds)) {
  sample_list <- c(sample_list, paste0(patient, "_", object_name[i]))
  assign(paste0(patient, "_", object_name[i]), readRDS(my_rds[i]))
}

# run get adjustedsegments and get output, name as $sample
# medicc2 parameter info will be stored
medicc2_params <- "medicc2_params.txt"
for (i in seq_along(object_name)) {
  # get ace output and best fits using minima
  current_object <- get(paste0(patient, "_", object_name[i]))
  model1 <- singlemodel(current_object, QDNAseqobjectsample = 1, ploidy = input_ploidys[i])
  # model1 is an object with different minima, trouble is there are
  # options either bestfit or lastminima we need to be able to adjust ploidy
  bestfit1 <- model1$minima[tail(which(model1$rerror == min(model1$rerror)), 1)]
  # I need to work on this
  log_txt <- paste0(patient, "_", object_name[i], " cellularity: ", bestfit1,
                    ", ploidy: ", model1$ploidy, "\n")
  cat(log_txt, file = medicc2_params, append = TRUE)
  assign(paste0(patient, "_", object_name[i], "_segments"),
         getadjustedsegments(current_object,
                             QDNAseqobjectsample = 1,
                             ploidy = input_ploidys[i],
                             cellularity = bestfit1))
}

# add sample_name column and name as $sample_segments
for (my_object in object_name) {
  current_object <- get(paste0(patient, "_", my_object, "_segments"))
  current_object <- explode_ranges(current_object, 1000000)
  all_bins <- get(paste0(patient, "_", my_object))
  all_bins  <- as.data.frame(all_bins@assayData$copynumber)
  all_bins <- all_bins %>%
    rownames_to_column(var = "position") %>%
    select(position) %>%
    separate(position, into = c("Chromosome", "Start", "End"), sep = "[:-]") %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End))

  current_object <- left_join(all_bins, current_object,
                              by = c("Chromosome", "Start", "End"))

  current_object$sample_id <- paste0(patient, "_", my_object)

  current_object <- current_object %>%
    dplyr::select("sample_id", "Chromosome", "Start", "End", "Copies") %>%
    rename(chrom = Chromosome, start = Start, end = End)
  assign(paste0(patient, "_", my_object, "_segments"),
         current_object, envir = .GlobalEnv)
}

dataframes <- mget(ls(pattern = "segments"))
out_df <- bind_rows(dataframes)
out_df <- pivot_wider(out_df, names_from = "sample_id", values_from = "Copies")
out_df <- out_df[complete.cases(out_df), ]
my_names <- paste0(patient, "_", object_name)

out_df <- out_df %>%
  group_by(chrom, across(my_names)) %>%
  summarize(
    start = min(start),
    end = max(end)
  )
out_df <- out_df %>%
  pivot_longer(
    cols = -c(chrom, start, end),    # Specify columns to gather
    names_to = "sample_id",      # Name for the new key column
    values_to = "Copies"    # Name for the new value column
  )
out_df$Diploid <- 2
out_df <- mutate(out_df, chrom = as.integer(chrom))
out_df <- out_df %>% arrange(chrom, start)

write.table(out_df, file = paste0(patient, ".tsv"), quote = FALSE,
            sep = "\t", row.names = FALSE)
