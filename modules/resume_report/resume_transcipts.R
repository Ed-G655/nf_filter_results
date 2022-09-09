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

#Calculate filtered total targets and transcripts
Total_venn_data.df <- Total_venn_data.df %>% 
  mutate(total_filtered_targets = total_input_targets - total_output_targets ) %>% 
  mutate(total_retained_targets_percent = (total_output_targets*100)/total_input_targets ) %>% 
  mutate(total_filtered_targets_percent = (total_filtered_targets*100)/total_input_targets ) %>% 
  mutate(total_filtered_transcripts = total_input_transcripts - total_output_transcripts ) %>% 
  mutate(total_retained_transcripts_percent = (total_output_transcripts*100)/total_input_transcripts ) %>% 
  mutate(total_filtered_transcripts_percent = (total_filtered_transcripts*100)/total_input_transcripts )

  
#Write dataframe of data by chromosome
write.table(venn_data.df, 
            file = "log_report.tsv", 
            sep = "\t", 
            row.names = F, 
            col.names = T)


#Write dataframe of data by chromosome
write.table(Total_venn_data.df, 
            file = "log_total_report.tsv", 
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
