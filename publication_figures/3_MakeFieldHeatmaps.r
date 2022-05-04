#! /usr/bin/Rscript

library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("~/Documents/Papers/1_Submitted/PUNNURI_Sorghum Aphid Resistance/GWAS_github/publication_figures/")
#setwd("C:/Users/jgwall/Desktop/FieldHeatmaps")
map2019=read.csv("FieldData/2019_map.csv", check.names=FALSE, colClasses="character")
map2020=read.csv("FieldData/2020_map.csv", check.names=FALSE, colClasses="character")
entries2019=read.csv("FieldData/2019_entry_list.csv", colClasses="character")
entries2020=read.csv("FieldData/2020_entry_list.csv", colClasses="character")
phenos=read.csv("FieldData/aphid damage and other Field data.csv", colClasses="character", check.names=FALSE)

# Convert map to long
names(map2019)[1]=names(map2020)[1]= "range"
layout2019 = pivot_longer(map2019, cols=names(map2019)[-1], names_to="row", values_to="plot")
layout2020 = pivot_longer(map2020, cols=names(map2020)[-1], names_to="row", values_to="plot")

# Merge into one data frame
## 2019
phenos$key=paste(phenos$PI, phenos$Rep, sep="_")
field2019 = entries2019 %>%
  mutate(key=paste(PI.Number, Rep, sep="_")) %>%
  inner_join(layout2019, by=c("X2019.Field"="plot")) %>%
  inner_join(phenos, by="key")
## 2020
field2020 = entries2020 %>%
  mutate(key=paste(PI, Rep, sep="_")) %>%
  inner_join(layout2020, by=c("X2020.Field" = "plot")) %>%
  inner_join(phenos, by="key")


# Helper function to make each plot identically
make_plot=function(myfield, pheno, name, legend_label, lowcolor, highcolor, limits=NULL){
  mydata= myfield[,c("row", "range", pheno)] %>%
    mutate_all(as.numeric)
  names(mydata)[3] = "pheno"
  
  myplot = ggplot(mydata) +
    aes(x=range, y=row, fill=pheno) +
    geom_raster() +
    ggtitle(name) +
    labs(x="Range", y="Row", fill=legend_label) +
    scale_fill_gradient(low=lowcolor, high=highcolor, limits=limits) +
    theme(axis.text=element_blank(), axis.ticks=element_blank(), 
          panel.background = element_blank(), title=element_text(face="bold"))
}

# Define colors; thanks to Paul Tol's notes on color schemes
blue_low="#C2E4EF"
blue_high="#364B9A"
red_low="#FEDA8B"
red_high="#A50026"
purp_low="#E7D4E8"
purp_high="#762A83"

# Get limits of count and damage scores
count.limits=c(field2019$AV.AC_Aug14_2019, field2019$AV.AC_Aug_28_2019,
               field2020$`AV.AC_Aug11-2020`, field2020$`AV.AC_Aug18-2020`) %>%
  as.numeric() %>%
  range(na.rm=T)
damage.limits=c(field2019$`2019D1`, field2019$`2019D2`, 
                field2020$`2020D1`, field2020$`2020D2`)%>%
  as.numeric() %>%
  range(na.rm=T)

# Make plots - 2019
flower.19=make_plot(field2019, "2019_Flowering", "Flowering Time (2019)", "DAP", blue_low, blue_high)
count1.19=make_plot(field2019, "AV.AC_Aug14_2019", "AC1 (2019)", "# Aphids", purp_low, purp_high, count.limits)
count2.19=make_plot(field2019, "AV.AC_Aug_28_2019", "AC2 (2019)", "# Aphids", purp_low, purp_high, count.limits)
d1.19=make_plot(field2019, "2019D1", "D1 (2019)", "Score", red_low, red_high, damage.limits)
d2.19=make_plot(field2019, "2019D2", "D2 (2019)", "Score", red_low, red_high, damage.limits)

# Export- 2019
png("Figure - Field Heatmaps 2019.png", width=8, height=6, units="in", res=300)
  blank = grid.rect(gp=gpar(col="white"))
  grid.arrange(flower.19, blank, count1.19, d1.19, count2.19, d2.19, ncol=2)
dev.off()

# Make plots - 2020
flower.20=make_plot(field2020, "2020_Flowering", "Flowering Time (2020)", "DAP", blue_low, blue_high)
count1.20=make_plot(field2020, "AV.AC_Aug11-2020", "AC1 (2020)", "# Aphids", purp_low, purp_high, count.limits)
count2.20=make_plot(field2020, "AV.AC_Aug18-2020", "AC2 (2020)", "# Aphids", purp_low, purp_high, count.limits)
d1.20=make_plot(field2020, "2020D1", "D1 (2020)", "Score", red_low, red_high, damage.limits)
d2.20=make_plot(field2020, "2020D2", "D2 (2020)", "Score", red_low, red_high, damage.limits)

# Export- 2020
png("Figure - Field Heatmaps 2020.png", width=8, height=6, units="in", res=300)
  blank = grid.rect(gp=gpar(col="white"))
  grid.arrange(flower.20, blank, count1.20, d1.20, count2.20, d2.20, ncol=2)
dev.off()