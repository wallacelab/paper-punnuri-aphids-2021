library(dplyr)
library(tidyr)


###########################################################3

#This function is used to find the candiate gene within 200kb region of a marker 
find_close_marker <- function(marker_list,gene_list,output_file_name){
  
  gene_list_row <- nrow(gene_list)
  
  for(i in 1:gene_list_row){
    ##select gene chromosome 
    chr_number = gene_list[i,"Chr"]
    ## select markers that are on the same chromosome with genes 
    temp_marker_list = subset(marker_list,Chr==chr_number)
    distance_min = Inf
    for(j in 1:nrow(temp_marker_list)){
      #select gene boundary
      start_position = gene_list[i,"Start"]
      stop_position = gene_list[i,"Stop"]
      #select marker position 
      marker_position = temp_marker_list[j,3]
      #calculate distance from marker position to the start and stop position of gene
      distance1 = abs(start_position - marker_position)
      distance2 = abs(stop_position - marker_position)
      ##calculate the range
        upper_bound=stop_position + 200000
        lower_bound=start_position - 200000
        
      ##determine the lower distance value
      if(marker_position >= lower_bound && marker_position <= upper_bound){
        if(distance1 < distance2){
          distance_now = distance1
        }
        else{
          distance_now = distance2
        }
        #record down the distance
        temp_marker_list$distance <- distance_now
        #bind the marker information and gene information
        best_hit = dplyr::bind_cols(temp_marker_list[j,],gene_list[i,])
        ##output the results
        write.table(best_hit,output_file_name,row.names = FALSE,col.names = FALSE,sep="\t",append = TRUE)
      }
    }
    }}

###########################################3


find_region <- function(candidate_gene_list){
  
  # make list of choromosome
  Chr_list <- unique(candidate_gene_list$Chr)
  
  #create an initial output list
  duplicate_marker_list <- candidate_gene_list[1,c(1:4)]
  
  #loop through chromosome
  for(i in Chr_list){
    chr_number = i
    #select genes on same chromosome
    temp_gene_list = subset(candidate_gene_list,Chr==chr_number)
    temp_gene_list = temp_gene_list[,c(1:4)]
    #select marker that repsent a region that contain more than 10 candidate genes 
    temp_gene_list <-unique(temp_gene_list %>% group_by(Marker_Name) %>% filter(n()>10))
    #add the result to output list
    duplicate_marker_list <- dplyr::bind_rows(duplicate_marker_list,temp_gene_list)
  }
  
  print(duplicate_marker_list)
  return(duplicate_marker_list)

}


###########################

########read in annotation file 
annotation <- read.csv("/home/hl46161/shorgun_project/data/formatted_Sbicolor_454_v3.1.1.defline.txt",sep="\t")
#marker_gene_annotation_list <- dplyr::inner_join(marker_gene_list,annotation,by="Name")

###########################
## read in global response data supplementary table 2
#global_response <- read.csv("/home/hl46161/shorgun_project/Genes/Table_4_Global Responses of Resistant and Susceptible Sorghum (Sorghum bicolor) to Sugarcane Aphid (Melanaphis sacchari).tsv",sep="\t")
global_response <- read.csv("/home/hl46161/shorgun_project/Genes/Table_2_Global Responses of Resistant and Susceptible Sorghum (Sorghum bicolor) to Sugarcane Aphid (Melanaphis sacchari).tsv",sep="\t")
## delete mid part
global_response <- global_response[-c(1:3),-c(6:29)]
## use first row to rename the column name
names(global_response) <- as.matrix(global_response[1, ])
## remove first row 
global_response <- global_response[-1, ]


#marker_gene_annotation_supplement <-  dplyr::inner_join(marker_gene_annotation_list,global_response,by = c("Name" = "Gene"))
#write.table(marker_gene_annotation_supplement,file = "/home/hl46161/shorgun_project/FarmCPU_without_flowering_time_marker_gene_annotation_overlap.txt", sep = "\t",row.names = FALSE)

################################################

## read in gene information of shorgun in gff format 
gff_file <- read.csv("/home/hl46161/shorgun_project/named_formatted_master_gene_list.gff3", sep = "\t")

## read in log2fold change information
#log2fold_information_global_response <- read.csv("/home/hl46161/shorgun_project/Genes/Table_1_Global Responses of Resistant and Susceptible Sorghum (Sorghum bicolor) to Sugarcane Aphid (Melanaphis sacchari).tsv",sep="\t",skip = 1)
#log2fold_information_global_response <- log2fold_information_global_response[,c("Gene","LFC")]

#gene annotation file 
annotation <- read.csv("/home/hl46161/shorgun_project/data/formatted_Sbicolor_454_v3.1.1.defline.txt",sep="\t")

## merge global response and gff file to get chr, start and stop condon information
merged_global_response <- dplyr::inner_join(gff_file,global_response,by = c("Name" = "Gene"))
#merged_global_response <- dplyr::inner_join(merged_global_response,log2fold_information_global_response,by = c("Name" = "Gene"))
merged_global_response <- dplyr::inner_join(merged_global_response,annotation,by = "Name")


#write.table(merged_global_response,file = "/home/hl46161/shorgun_project/formated_DE_gene.txt", sep = "\t",row.names = FALSE)

######################################################


#read in BMC datasets 
BMC_genomics_2_wk_resistant <- read.csv("/home/hl46161/shorgun_project/Genes/12864_2018_5095_MOESM2_ESM.tsv",sep="\t",skip=6)
BMC_genomics_2_wk_susceptible  <- read.csv("/home/hl46161/shorgun_project/Genes/12864_2018_5095_MOESM3_ESM.tsv",sep="\t",skip=6)
BMC_genomics_6_wk_resistant  <- read.csv("/home/hl46161/shorgun_project/Genes/12864_2018_5095_MOESM4_ESM.tsv",sep="\t",skip=5)
BMC_genomics_6_wk_susceptible  <- read.csv("/home/hl46161/shorgun_project/Genes/12864_2018_5095_MOESM5_ESM.tsv",sep="\t",skip=5)

## delete mid part
colnames(BMC_genomics_2_wk_resistant)
BMC_genomics_2_wk_resistant <- BMC_genomics_2_wk_resistant[c(1:(nrow(BMC_genomics_2_wk_resistant)-4)),c("Feature.ID","Descriptiona","Fold.changeb","FDR.p.value","Categoryc")]
BMC_genomics_2_wk_susceptible <- BMC_genomics_2_wk_susceptible[c(1:(nrow(BMC_genomics_2_wk_susceptible)-4)),c("Feature.ID","Descriptiona","Fold.changeb","FDR.p.value","Categoryc")]
BMC_genomics_6_wk_resistant <- BMC_genomics_6_wk_resistant[c(1:(nrow(BMC_genomics_6_wk_resistant)-4)),c("Feature.ID","Descriptiona","Fold.changeb","FDR.p.value","Categoryc")]
BMC_genomics_6_wk_susceptible <- BMC_genomics_6_wk_susceptible[c(1:(nrow(BMC_genomics_6_wk_susceptible)-4)),c("Feature.ID","Descriptiona","Fold.changeb","FDR.p.value","Categoryc")]

#rename the BMC genomics log fold change column tp specify experiment time 
colnames(BMC_genomics_2_wk_resistant)[3] <- "Foldchange_resistant_2wk"

colnames(BMC_genomics_2_wk_susceptible)[3] <- "Foldchange_susceptible_2wk"

colnames(BMC_genomics_6_wk_resistant)[3] <- "Foldchange_resistant_6wk"

colnames(BMC_genomics_6_wk_susceptible)[3] <- "Foldchange_susceptible_6wk"


###
#join the seperate BMC dataset tpgether
BMC_genomics <- dplyr::full_join(BMC_genomics_2_wk_resistant,BMC_genomics_2_wk_susceptible,by="Feature.ID",copy=TRUE)
BMC_genomics <- BMC_genomics[,c("Feature.ID","Descriptiona.x","Foldchange_susceptible_2wk","Foldchange_resistant_2wk","FDR.p.value.x","Categoryc.x")]

BMC_genomics <- dplyr::full_join(BMC_genomics,BMC_genomics_6_wk_susceptible,by="Feature.ID",copy=TRUE)
BMC_genomics <- BMC_genomics[,c("Feature.ID","Descriptiona.x","Foldchange_susceptible_2wk","Foldchange_resistant_2wk","Foldchange_susceptible_6wk","FDR.p.value.x","Categoryc.x")]

BMC_genomics <- dplyr::full_join(BMC_genomics,BMC_genomics_6_wk_resistant,by="Feature.ID",copy=TRUE)
BMC_genomics <- BMC_genomics[,c("Feature.ID","Descriptiona.x","Foldchange_susceptible_2wk","Foldchange_resistant_2wk","Foldchange_susceptible_6wk","Foldchange_resistant_6wk","FDR.p.value.x","Categoryc.x")]


#merge the BMC datasets to get chr, start and stop information 
merged_BMC_genomics <- dplyr::inner_join(gff_file,BMC_genomics,by = c("Name" = "Feature.ID"))
merged_BMC_genomics <- dplyr::inner_join(merged_BMC_genomics,annotation,by = "Name")

####################################


########
#read in without flowering time marker list 
final_marker_list <- read.csv("/home/hl46161/shorgun_project/final_hit/combined.unique_final_list.tsv",sep="\t",col.names = c("Marker_Name","Chr","Position",'traits',"Source"))

########################

#
find_close_marker(final_marker_list,merged_global_response,"/home/hl46161/shorgun_project/final_hit/final_hit_200k_range_marker_gene_list_global_response.txt")

## read in to give column name 
final_hit_200k_range_marker_gene_list_global_response <- read.csv("/home/hl46161/shorgun_project/final_hit/final_hit_200k_range_marker_gene_list_global_response.txt",sep="\t",
                                                                  col.names = c("Marker_Name","Chr",'Position',"Traits","Source_dataset","Distance","Chr","Source","Type","Start","Stop","Score","Strand","Phase","ID","Name","AncestorIdentifier","Genotype comparison",
                                                                                "Day","Expression","Number of genes in VD","Best hit Arabidopsis gene","Arabidopsis definition","Best hit Oryza sativa","Oryza sativa definition","Reference_source","function"))


write.table(final_hit_200k_range_marker_gene_list_global_response,"/home/hl46161/shorgun_project/final_hit/final_hit_200k_range_marker_gene_list_global_response.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)




########## find marker in 200 K range 

find_close_marker(final_marker_list,merged_BMC_genomics,"/home/hl46161/shorgun_project/final_hit/final_hit_200k_range_marker_gene_list_BMC_genomics.txt")

##################

final_marker_list_200k_range_marker_gene_list_BMC_genomics <- read.csv("/home/hl46161/shorgun_project/final_hit/final_hit_200k_range_marker_gene_list_BMC_genomics.txt",sep="\t",
                                                                                   col.names = c("Marker_Name","Chr","Position",'traits',"Source_dataset","Distance","Chr","Source","Type","Start","Stop","Score","Strand","Phase","ID","Name","AncestorIdentifier",
                                                                                                 "Descriptiona.x","Foldchange_susceptible_2wk","Foldchange_resistant_2wk","Foldchange_susceptible_6wk","Foldchange_resistant_6wk","FDR.p.value.x","Categoryc.x","Reference_source","function"))

write.table(final_marker_list_200k_range_marker_gene_list_BMC_genomics,"/home/hl46161/shorgun_project/final_hit/final_hit_200k_range_marker_gene_list_BMC_genomics.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)


##filter out the marker information
unique_final_marker_list_200k_range_marker_gene_list_BMC_genomics <- final_marker_list_200k_range_marker_gene_list_BMC_genomics[,c(7:ncol(final_marker_list_200k_range_marker_gene_list_BMC_genomics))]
##remove diplicate results
unique_final_marker_list_200k_range_marker_gene_list_BMC_genomics <- unique(unique_final_marker_list_200k_range_marker_gene_list_BMC_genomics)
unique_final_marker_list_200k_range_marker_gene_list_BMC_genomics

##filter out the marker information
unique_final_hit_200k_range_marker_gene_list_global_response <- final_hit_200k_range_marker_gene_list_global_response[,c(7:ncol(final_hit_200k_range_marker_gene_list_global_response))]
##remove diplicate results
unique_final_hit_200k_range_marker_gene_list_global_response <- unique(unique_final_hit_200k_range_marker_gene_list_global_response)
unique_final_hit_200k_range_marker_gene_list_global_response

#output the results
write.table(unique_final_marker_list_200k_range_marker_gene_list_BMC_genomics,"/home/hl46161/shorgun_project/final_hit/unique_final_hit_200k_range_marker_gene_list_BMC_genomics.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)
write.table(unique_final_hit_200k_range_marker_gene_list_global_response,"/home/hl46161/shorgun_project/final_hit/unique_final_hit_200k_range_marker_gene_list_global_response.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)

unique(final_marker_list_200k_range_marker_gene_list_BMC_genomics$Name)
unique(final_hit_200k_range_marker_gene_list_global_response$Name)

####################################

##find region with more than 10 candiate genes 
global_response_duplicate_region <- find_region(final_hit_200k_range_marker_gene_list_global_response)
BMC_genomics_region <- find_region(final_marker_list_200k_range_marker_gene_list_BMC_genomics)

write.table(global_response_duplicate_region,"/home/hl46161/shorgun_project/final_hit/global_response_region.txt",row.names = FALSE,col.names = TRUE,sep="\t")
write.table(BMC_genomics_region,"/home/hl46161/shorgun_project/final_hit/BMC_genomics_region.txt",row.names = FALSE,col.names = TRUE,sep="\t")

###############################################################
