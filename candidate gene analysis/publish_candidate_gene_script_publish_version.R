library(dplyr)
library(tidyr)


###########################################################3

#This function is used to find the candiate gene within 200kb region of a marker 
find_close_marker <- function(marker_list,gene_list){
  
  gene_list_row <- nrow(gene_list)
  colnames(marker_list)[2] = "Marker_Chr"
  marker_list$distance <- 0
  best_hit = dplyr::bind_cols(marker_list[1,],gene_list[1,])
  #print(best_hit)
  for(i in 1:gene_list_row){
    ##select gene chromosome 
    chr_number = gene_list[i,"Chr"]
    ## select markers that are on the same chromosome with genes 
    temp_marker_list = subset(marker_list,Marker_Chr==chr_number)
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
        temp_best_hit = dplyr::bind_cols(temp_marker_list[j,],gene_list[i,])
        best_hit <- dplyr::bind_rows(best_hit,temp_best_hit)

      }
    }
  }
  best_hit <- best_hit[-1,]
  return(best_hit)
}

#final_hit_200k_range_marker_gene_list_global_response_full <- find_close_marker(final_marker_list,merged_global_response)


###########################################3


find_region <- function(candidate_gene_list){
  
  # make list of choromosome
  Chr_list <- unique(candidate_gene_list$Chr)
  
  #create an initial output list
  duplicate_marker_list <- candidate_gene_list[1,c("Marker_Name","Chr","Position","Name")]
  
  #loop through chromosome
  for(i in Chr_list){
    chr_number = i
    #select genes on same chromosome
    temp_gene_list = subset(candidate_gene_list,Chr==chr_number)
    temp_gene_list = unique(temp_gene_list[,c("Marker_Name","Chr","Position","Name")])
  
    #select marker that repsent a region that contain more than 10 candidate genes 
    temp_gene_count_list <- unique(temp_gene_list %>% group_by(Marker_Name) %>% filter(n()>10) %>% summarize(Chr=Chr,Position=Position,candidate_gene_count=n()))
    print(temp_gene_count_list)
    duplicate_marker_list <- dplyr::bind_rows(duplicate_marker_list,temp_gene_count_list)
  }
  duplicate_marker_list <- duplicate_marker_list[-1,]
  #print(duplicate_marker_list)
  return(duplicate_marker_list)
}

#test_final_hit_200k_range_marker_gene_list_global_response <- find_close_marker(final_marker_list,merged_BMC_genomics)

###########################

################################################

## read in gene information of shorgun in gff format 
gff_file <- read.csv("/home/hl46161/shorgun_project/named_formatted_master_gene_list.gff3", sep = "\t")

#gene annotation file 
annotation <- read.csv("/home/hl46161/shorgun_project/data/formatted_Sbicolor_454_v3.1.1.defline.txt",sep="\t")

global_response <- read.csv("/home/hl46161/shorgun_project/Genes/Table_3_Global Responses of Resistant and Susceptible Sorghum (Sorghum bicolor) to Sugarcane Aphid (Melanaphis sacchari).csv",sep="\t",skip=1)

#filter based on module 
global_response <- subset(global_response,Module.Assignment %in% c("1","12","18","23","4","9","10","14","16","17"))

global_response <- global_response[-c(3:8)]

## merge global response and gff file to get chr, start and stop condon information
merged_global_response <- dplyr::inner_join(gff_file,global_response,by = c("Name" = "Gene"))
#merged_global_response <- dplyr::inner_join(merged_global_response,log2fold_information_global_response,by = c("Name" = "Gene"))
merged_global_response <- dplyr::inner_join(merged_global_response,annotation,by = "Name")
merged_global_response <- subset(merged_global_response,Chr %in% c(1:10))

write.table(merged_global_response,file = "/home/hl46161/shorgun_project/final_hit/global_response_filtered_DEG.tsv", sep = "\t",row.names = FALSE)

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
merged_BMC_genomics <- subset(merged_BMC_genomics,Chr %in% c(1:10))

write.table(merged_BMC_genomics,"/home/hl46161/shorgun_project/final_hit/BMC_DEG_list.tsv",row.names = FALSE,col.names = TRUE,sep="\t")

####################################


########
#read in without flowering time marker list 
final_marker_list <- read.csv("/home/hl46161/shorgun_project/final_hit/combined.unique_final_list.tsv",sep="\t",col.names = c("Marker_Name","Chr","Position",'traits',"Source"))

########################

#run the candidate gene search for global response in 200k range 
final_hit_200k_range_marker_gene_list_global_response <- find_close_marker(final_marker_list,merged_global_response)
final_hit_200k_range_marker_gene_list_global_response 
#output the results 
write.table(final_hit_200k_range_marker_gene_list_global_response,"/home/hl46161/shorgun_project/final_hit/final_hit_200k_range_marker_gene_list_global_response.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)

##########
#run the candidate gene search for BMC genomics in 200k range 
final_hit_200k_range_marker_gene_list_BMC_genomics <- find_close_marker(final_marker_list,merged_BMC_genomics)
final_hit_200k_range_marker_gene_list_BMC_genomics

write.table(final_hit_200k_range_marker_gene_list_BMC_genomics,"/home/hl46161/shorgun_project/final_hit/final_hit_200k_range_marker_gene_list_BMC_genomics.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)


##filter out the marker information
unique_final_hit_200k_range_marker_gene_list_BMC_genomics <- final_hit_200k_range_marker_gene_list_BMC_genomics[,c(7:ncol(final_hit_200k_range_marker_gene_list_BMC_genomics))]
##remove diplicate results
unique_final_hit_200k_range_marker_gene_list_BMC_genomics <- unique(unique_final_hit_200k_range_marker_gene_list_BMC_genomics)
unique_final_hit_200k_range_marker_gene_list_BMC_genomics

##filter out the marker information
unique_final_hit_200k_range_marker_gene_list_global_response <- final_hit_200k_range_marker_gene_list_global_response[,c(7:ncol(final_hit_200k_range_marker_gene_list_global_response))]
##remove diplicate results
unique_final_hit_200k_range_marker_gene_list_global_response <- unique(unique_final_hit_200k_range_marker_gene_list_global_response)
unique_final_hit_200k_range_marker_gene_list_global_response


#output the results
write.table(unique_final_hit_200k_range_marker_gene_list_BMC_genomics,"/home/hl46161/shorgun_project/final_hit/unique_final_hit_200k_range_marker_gene_list_BMC_genomics.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)
write.table(unique_final_hit_200k_range_marker_gene_list_global_response,"/home/hl46161/shorgun_project/final_hit/unique_final_hit_200k_range_marker_gene_list_global_response.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)


####################################

##find region with more than 10 candiate genes 
global_response_region <- find_region(final_hit_200k_range_marker_gene_list_global_response)
BMC_genomics_region <- find_region(final_hit_200k_range_marker_gene_list_BMC_genomics)

write.table(global_response_region,"/home/hl46161/shorgun_project/final_hit/global_response_region.txt",row.names = FALSE,col.names = TRUE,sep="\t")
write.table(BMC_genomics_region,"/home/hl46161/shorgun_project/final_hit/BMC_genomics_region.txt",row.names = FALSE,col.names = TRUE,sep="\t")

###############################################################


## read in global response data supplementary table 2
global_response <- read.csv("/home/hl46161/shorgun_project/Genes/Table_3_Global Responses of Resistant and Susceptible Sorghum (Sorghum bicolor) to Sugarcane Aphid (Melanaphis sacchari).csv",sep="\t",skip=1)

global_response_full <- global_response[-c(3:8)]

#merged with gff file to get location information
full_merged_global_response <- dplyr::inner_join(gff_file,global_response_full,by = c("Name" = "Gene"))
## add annotation 
full_merged_global_response <- dplyr::inner_join(full_merged_global_response,annotation,by = "Name")
#filter out chromosome that are not on Chr 1 to 10 
full_merged_global_response <- subset(full_merged_global_response,Chr %in% c(1:10))
write.table(full_merged_global_response,file = "/home/hl46161/shorgun_project/final_hit/global_response_full_DEG.tsv", sep = "\t",row.names = FALSE)


full_final_hit_200k_range_marker_gene_list_global_response <- find_close_marker(final_marker_list,full_merged_global_response)
full_final_hit_200k_range_marker_gene_list_global_response 
#output the results 
write.table(full_final_hit_200k_range_marker_gene_list_global_response,"/home/hl46161/shorgun_project/final_hit/full_final_hit_200k_range_marker_gene_list_global_response.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)

##filter out the marker information
unique_final_hit_200k_range_marker_gene_list_global_response_full <- full_final_hit_200k_range_marker_gene_list_global_response[,c(7:ncol(full_final_hit_200k_range_marker_gene_list_global_response))]
##remove diplicate results
unique_final_hit_200k_range_marker_gene_list_global_response_full <- unique(unique_final_hit_200k_range_marker_gene_list_global_response_full)
unique_final_hit_200k_range_marker_gene_list_global_response_full

write.table(unique_final_hit_200k_range_marker_gene_list_global_response_full,"/home/hl46161/shorgun_project/final_hit/full_unique_final_hit_200k_range_marker_gene_list_global_response_full.txt",row.names = FALSE,sep="\t",na = "",col.names = TRUE)

#find high genes count (>10) region 
full_global_response_duplicate_region <- find_region(full_final_hit_200k_range_marker_gene_list_global_response)
write.table(full_global_response_duplicate_region,"/home/hl46161/shorgun_project/final_hit/full_global_response_region.txt",row.names = FALSE,col.names = TRUE,sep="\t")

