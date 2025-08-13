setwd("/groups/berger/user/vikas.shukla/States_AT/Analysis/Landscapes/AP_clustering/")
library(RColorBrewer)
library(valr)
library(tidyverse)
library(dtwclust)
library(apcluster)
library(dendextend)

#-------------- Function to create a directory if it doesn't exist -------------
create_dir_if_needed <- function(dir_name) {
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
}
# Define and create the necessary directories
dir_rdata <- "Rdata"
dir_plots <- "Cluster_plots"
dir_subclusters <- "Merged_subclusters"
dir_bed <- "BED_files"
# Use the function to create directories
create_dir_if_needed(dir_rdata)
create_dir_if_needed(dir_plots)
create_dir_if_needed(dir_subclusters)
create_dir_if_needed(dir_bed)
#--------
# width of chromHMM state (200bp)
state_size=200 
## data files
gff_file = "data/TAIR10_GFF3_genes_transposons.gff"
state_bed= "data/AT_26_segments.bed"
emission_matrix= "data/emissions_26.txt"

#Defining the new names and colors for the states
cs <- c(brewer.pal(8,"Reds")[8:2],
            brewer.pal(9,"Blues")[c(6,8)],
            brewer.pal(7,"Purples")[c(3,5,7)],
            brewer.pal(7,"Greys")[6:3],
            brewer.pal(3,"YlOrBr")[2:1],
            brewer.pal(9,"Greens")[1:8])
map_to_paper <- data.frame(
  from = c("21","22","19","20","18","23","17",
               "12","11","10","8","9",
               "24","25","26","16",
               "6","7","5","4","3","2","15","1","14","13"),
  to = c("H1","H2","H3","H4","H5","H6","H7",
             "F1","F2","F3","F4","F5",
             "I1","I2","I3","I4",
             "E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"))
names(cs) <- map_to_paper$to

## Read in files, manipulate to right format
gff <- read_delim(gff_file,delim="\t",comment = "#",col_names = FALSE)
emiss_mat <- read_delim("data/emissions_26.txt")%>%dplyr::rename(state=`State (Emission order)`)
states <- read_delim("data/AT_26_segments.bed",col_names = c("chrom","start","end","state")) %>% 
  mutate(chrom=as.character(chrom))

## map to better state names
states=states%>%mutate(state=gsub("^E","",state))%>% 
  full_join(map_to_paper,by=c("state"="from"))%>%
  select(-state)%>%
  dplyr::rename(state=to)

emiss_mat = emiss_mat%>%mutate(state=as.character(state))%>% 
  full_join(map_to_paper,by=c("state"="from"))%>%
  select(-state)%>%
  dplyr::rename(state=to)%>%
  relocate(state) 

## Extract protein_coding_gene from gff 
gene_pos=gff%>%filter(!X1%in%c("ChrC","ChrM"))%>% # remove chloroplasm, mitocondria
  filter(grepl('protein_coding_gene', X9))%>% 
  mutate(geneId=str_match(X9,"^ID=([A-Za-z0-9]+);")[,2]) %>% 
  select(X1,X4,X5,X7,geneId)%>%
  dplyr::rename(chrom=X1,start=X4,end=X5,strand=X7)%>%
  mutate(chrom=str_match(chrom,"^Chr(.+)")[,2])

# clean up 
rm("gff","map_to_paper","emission_matrix","gff_file","state_bed")
