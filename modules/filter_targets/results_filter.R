#Load libraries
library("dplyr")
library("stringr")
library("vroom")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only

#args[1] <-"10"

#args[2] <-"Emsembl_Canonical.tsv"

#args[3] <- "chr10.mirmapout"

# GET CHROM
chromosome <- args[1]

# Load Canonical transcript IDs (ensemble canonical) database from args
Ensemble_Canonical <- args[2]

## get miRNA targets from args
miRNOme_out <- args[3]



# Load Canonical transcript IDs (ensemble canonical) database
Ensemble_Canonical.df <- vroom(Ensemble_Canonical)  %>% na.omit()

#Load miRNA targets
miRNOme_out.df <- vroom(miRNOme_out)

 #Rename first Column 
colnames(miRNOme_out.df)[1] <- "GeneID"

#Get transcript stable ID version
miRNOme_out.df <- miRNOme_out.df %>%  mutate( `Transcript stable ID` = str_match_all(GeneID, pattern = "ENST00000\\d+" ))

miRNOme_out.df <- miRNOme_out.df %>% mutate(`Transcript stable ID` = 
                                        as.character(`Transcript stable ID`))

#Get only canonical targets 
semi_join_miRNOme.df <-  semi_join(miRNOme_out.df, Ensemble_Canonical.df, 
                                by =  "Transcript stable ID")
                         
Output_file <- semi_join_miRNOme.df %>% select(-`Transcript stable ID`)


#Get file data

file_extencion <- str_match(miRNOme_out, pattern = "[^chr\\d+\\.].+")


#Write filtered targets
write.table(x = Output_file, 
            file = str_interp("${chromosome}.filtered.${file_extencion}"), 
            sep = "\t", 
            row.names =  F, 
            col.names = T)

### Write log file
input_targets  <- as.integer(nrow(miRNOme_out.df))
output_targets  <- nrow(semi_join_miRNOme.df)
filtered_targets <- input_targets - output_targets
retained_percent  <- (output_targets*100)/input_targets
filtered_percent <- (filtered_targets*100)/input_targets


### Get unique Transcripts IDs resume
Canonical_transcripts <- Ensemble_Canonical.df %>%  pull(`Transcript stable ID`) %>%  unique() %>% length()
input_transcripts <- miRNOme_out.df %>%  pull(`Transcript stable ID`) %>% unique() %>% length()
output_transcripts <- semi_join_miRNOme.df %>%  pull(`Transcript stable ID`) %>% unique() %>% length()
filtered_transcripts <- input_transcripts - output_transcripts
retained_transcripts_percent <- (output_transcripts*100)/input_transcripts
filtered_transcripts_percent <- (filtered_transcripts*100)/input_transcripts


#### Write log dataframe
log_file.df <- data_frame(chromosome,
                       input_targets, 
                       output_targets, 
                       filtered_targets, 
                       retained_percent, 
                       filtered_percent,
                       input_transcripts,
                       output_transcripts,
                       filtered_transcripts,
                       retained_transcripts_percent,
                       filtered_transcripts_percent)


#Write filtered targets
write.table(x = log_file.df, 
            file = str_interp("${miRNOme_out}.log"), 
            sep = "\t", 
            row.names =  F, 
            col.names = T)


