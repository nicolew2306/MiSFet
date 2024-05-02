### Loading the libraries
##### Need devtools to be installed through conda. This script should be run from within that conda env. Checking if they're already installed, if not install them. Note these need to be installed within a devtools conda env.


if (!require("devtools", character.only = TRUE)) {  
  warning("You may need to use a custom conda environment. Please look at the readme file.")    
} else {
  library(FEAST)
}


package_list <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes","devtools")

for (package_name in package_list) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name)
    library(package_name, character.only = TRUE)
  } else {
    library(package_name, character.only = TRUE)
  }
}

if (!require("FEAST", character.only = TRUE)) {
  # need to install with devtools from github
  devtools::install_github("cozygene/FEAST")
  library(FEAST)
} else {
  library(FEAST)
}

#### signature output data (directory) and signature metadata (file) paths input by user. 

if (length(commandArgs(trailingOnly = TRUE)) < 2) {
  stop("Provide paths for both the data directory and metadata file. I can't read your mind.")
}

# Get the paths from the command-line arguments
data_path <- commandArgs(trailingOnly = TRUE)[1]
metadata_path <- commandArgs(trailingOnly = TRUE)[2]


#### Reading in the sink_source metadata file from the path given above.

### possibility to enter the paths directly into the script
#metadata_path <- "/path/to/sink_source.csv"
#data_path <- "/PATH/to/species"

## this is the same metadata file that is used for Signature. Make sure the format matches that exactly.
sink_source_info <- read.csv(metadata_path,sep=",")

##checking to make sure there is no "_" in the family_id that would complicate the matching to samples.
if (any(grepl("_", sink_source_info$family_id))) {
    stop("Please rename you family_id column to not include underscores.")
}

#### Here we build the FEAST metadata file from the sink_source Signature file.

## reorganize the sink_source so that the FEAST metadata can be built from it.
rownames(sink_source_info) <- sink_source_info$family_id
sink_source_info <- sink_source_info[,-c(1)]

## transposing for building the dataframe
sink_source_info_t <- t(sink_source_info)
sink_source_info_t <- data.frame(sink_source_info_t)

## prepping the column dataframes to build the FEAST metadata file

feast_metadata_sampleid <- feast_metadata_env <- feast_metadata_id <- feat_metadata_sinksource <- NULL

for (i in 1:ncol(sink_source_info_t)) {
    
    ## SampleID column
    feast_metadata_sampleid <- rbind(feast_metadata_sampleid,data.frame(sink_source_info_t[,i]))
    ## Env column
    feast_metadata_env <- rbind(feast_metadata_env,data.frame(rownames(sink_source_info_t)))
    ## sinksouce columns
    SourceSink_col <- data.frame(c("sink", rep("source", times = nrow(sink_source_info_t) - 1)))
    feat_metadata_sinksource <- rbind(feat_metadata_sinksource,SourceSink_col)
    ## id column  
    id_col <- data.frame(rep(colnames(sink_source_info_t)[i],times=nrow(sink_source_info_t)))
    feast_metadata_id <- rbind(feast_metadata_id,id_col)
}  

## combining and renaming the columns correctly
feast_metadata <- cbind(feast_metadata_sampleid,feast_metadata_env,feat_metadata_sinksource,feast_metadata_id)
colnames(feast_metadata) <- c("SampleID","Env","SourceSink","id")

#### in case the names contain "-" since R replaces them with ".". If there are no "-", then it will do nothing.
feast_metadata$Env <- sub("\\.", "-", feast_metadata$Env)

#### individual_sink_source_feast_prep function for building the required dataframes. It creates two list of dataframes each containg the otu_count and the metadata along with the sample name for a specific sink and its sources. 



individual_sink_feast_prep <- function(file) {
  # Extracting the file name from the who path
  parts_file_name <- strsplit(file,"/")
  file_name <- parts_file_name[[1]][length(parts_file_name[[1]])]  
    
  # This is the Signature output so it should be fine. Unless they name their samples with "_"
  sample_parts <- strsplit(file_name, "_") 
  # The names have to match the metadata file, so as long as they don't use "_" in their sample_ids it should
  ## be ok.
  sample_name <- sample_parts[[1]][length(sample_parts[[1]])-1] 
    
  # Read the CSV file into a data frame
  df <- read.csv(file, sep = ",", header = FALSE, row.names = 1)
  df[is.na(df)] <- 0
  
  # Assigning column names
  col_names <- sink_source_info[rownames(sink_source_info) == sample_name, ]
  colnames(df) <- col_names
  
  # Removing rows with all 0s to avoid a "0 depth" error
  df_t <- t(df)
  df_t <- df_t[rowSums(df_t != 0, na.rm = TRUE) > 0, ]
  
  metadata_df <- feast_metadata[feast_metadata$id == sample_name, ]
  rownames(metadata_df) <- metadata_df$SampleID
  metadata_df <- metadata_df[, -which(names(metadata_df) == "SampleID")]

  
  # Return a list containing the sample name, data frame, and metadata
  return(list(sample_name = sample_name, data_frame = df_t, metadata = metadata_df))
}


#### output from Signature for individual prep run. Converting to input using the individual_sink_feast_prep function.

#### this is for individual source only.

file_list <- list.files(path=data_path, pattern="counts.csv$", all.files=T, full.names=T)

# lists to store results
count_list <- list()
metadata_list <- list()


for (file in file_list) {
  result <- individual_sink_feast_prep(file)
  count_list[[result$sample_name]] <- result$data_frame
  metadata_list[[result$sample_name]] <- result$metadata
}


#### Starting the FEAST run.
result_path <- "./"
#### looping from input file
print("This might take some time. Be patient. It's good for the soul...occasionally.")
if (length(count_list) == length(metadata_list)) {
  
  # Iterate through the lists using a loop
  for (i in seq_along(count_list)) {
    print(i)
    otus <- count_list[[i]]
    metadata <- metadata_list[[i]]
    ## running FEAST
    FEAST(C = otus, metadata = metadata , different_sources_flag = 1, dir_path = result_path,
                      outfile=metadata_list[[i]][1,3])    
  }
    
} else {
  warning("The lists have different lengths.")

}

metadata_file_path <- paste(metadata_path,"feast_metadata.txt",sep="_")
write.table(metadata,metadata_file_path)

print("Done!")
