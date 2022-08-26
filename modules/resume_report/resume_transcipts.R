## Cargamos pacman
library( "pacman" )

#cargamos mas paquetes
p_load( "vroom",
        "dplyr",
        "purrr")




# Cargamos todos los tsv en un solo dataframe
filenames <- list.files( pattern="*.log" )

venn_data.df <- map_df( filenames, vroom)

names(venn_data.df)

#Sum colums
Total_venn_data.df <- venn_data.df %>%  summarise( total_input_targets= sum(input_targets)) %>% 
  bind_cols(venn_data.df %>% summarise( total_output_targets = sum(output_targets))) %>% 
  bind_cols(venn_data.df %>%  summarise( total_input_transcripts = sum(input_transcripts))) %>% 
  bind_cols(venn_data.df %>%  summarise( total_output_transcripts = sum(output_transcripts))) 
  
#Write dataframe of data by chromosome
write.table(venn_data.df, 
            file = "changes_by_chr_genes.tsv", 
            sep = "\t", 
            row.names = F, 
            col.names = T)


# create pairwise Venn diagram
#resume_targets.p <- draw.pairwise.venn(area1=Total_venn_data.df$total_output_targets,
#                   area2=Total_venn_data.df$total_input_targets,
#                   cross.area=Total_venn_data.df$total_output_targets,
#                   category=c("Output Targets","Input Targets"),fill=c("Red","Yellow"))


# create pairwise Venn diagram
#draw.pairwise.venn(area1=Total_venn_data.df$total_input_transcripts,
 #                  area2=Total_venn_data.df$total_output_transcripts,
  #                 cross.area=Total_venn_data.df$total_output_transcripts,
   #                category=c("Input transcripts","Output transcripts"),fill=c("Red","Yellow"))
