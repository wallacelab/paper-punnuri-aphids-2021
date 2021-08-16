#! /usr/bin/Rscript

# Helper script including functions for making composite figures

require(ggplot2)
require(ggthemes)
require(scales)

# Load GLM data
load_glm = function(infile, debug=FALSE){
    if(!is.null(infile) && file.exists(infile)){
        cat("\tLoading GLM results from",infile,"\n")
        glm = read.csv(infile)
    }else {
        cat("\tUnable to locate GLM file",infile,"; creating dummy data\n")
        glm = data.frame(Chr=1, Pos=1, p=c(1, 0.1), perm_p=1)
    }
    glm$log.p = -log10(glm$p)  # -log10 of p-value
    
    # For debugging, take random sample of GLM results to speed process
    if(debug){
        set.seed(1);
        glm = glm[sample(1:nrow(glm), size=1000, replace=F),]
    }
    
    return(glm)
}

# Load FarmCPU RMIP data
load_farm = function(infile){
    if(file.exists(infile)){
        cat("\tLoading FARM-CPU results from",infile,"\n")
        farm = read.csv(infile)
    }else{
        cat("\tUnable to locate FarmCPU file",infile,"; creating dummy data\n")
        farm = data.frame(snp=character(), chr=integer(), pos=integer(), rmip=numeric(), significant=character())
    }
    names(farm)[names(farm)=="chr"] = "Chr"
    names(farm)[names(farm)=="pos"] = "Pos"
    return(farm)
}

# Load and format offsets
load_offsets = function(infile){
    offsets=read.csv(infile, row.names='chr')
    offsets$offset = offsets$offset + as.numeric(rownames(offsets)) * 5e6 # Add padding between chromosomes
    offsets$midpoint = offsets$offset + (offsets$length/2) # Calculate midpoint of each chromosome
    return(offsets)   
}

# Helper function to get x position for plotting a point given a vector of chrom and positions
get_offset = function(chr, pos){
    myoffset = offsets[as.character(chr),"offset"]
    return(myoffset + pos)
}


# Plot a single manhattan plot
plot_manhattan = function(trait_name, glm_file, farm_file, args){

    # Load in data
    cat("Plotting data for",trait_name, "\n")
    glm = load_glm(glm_file, debug=args$debug)
    farm = load_farm(farm_file)
    
    # Subset out permutation significance
    glm.sig = subset(glm, glm$perm_p <= args$perm_cutoff)
    
    # Reformat RMIP data
    farm.all = farm # Backup of complete RMIP data
    farm = subset(farm, farm$rmip >= args$rmip_cutoff)
    
    
    # Set up plotting colors
    chrom_colors=RColorBrewer::brewer.pal(3, "Set1")[1:2]
    my_alpha = 0.5
    if(all(glm$Pos==1)){ # Set alpha to 0 for GLM if is a dummy dataset (all positions at 1)
        my_alpha=0
    }
 
    # Basic plotting
    myplot = ggplot(glm, mapping=aes(x=get_offset(Chr, Pos), y=log.p, color=(Chr%%2 == 0))) + 
        geom_point(alpha=my_alpha) +  # Make a scatterplot
        scale_color_manual(values=chrom_colors) +
        theme_few() + # Low-clutter theme
        scale_x_continuous(breaks=offsets$midpoint, minor_breaks = offsets$offset, labels=rownames(offsets)) +  # Set labels and line positions for X axis
        theme(legend.position = "none") +  # Turn off legend
        #labs(x="Chromosome", y=expression(bold(paste("-", log[10], " p-value (GLM)")))) + # X and Y axis labels
        labs(x=element_blank(), y=element_blank()) + # X and Y axis labels
        ggtitle(trait_name) +
        theme(title = element_text(size=10, face='bold')) # Set title 
        
    
    # Add significant GLM hits
    myplot = myplot +
        geom_point(data=glm.sig, size=2, color='black', fill=alpha('white', 0), shape=21, stroke=1.5)
    
    # Adjust Y axis - add right axis and set significant digits for everything
    if(nrow(farm.all)>0){
        farm.max = max(farm.all$rmip)
    }else{
        farm.max=1
    }
    scale_trans = farm.max / max(glm$log.p)    
    myplot = myplot +
        scale_y_continuous(labels = scales::number_format(accuracy=0.1),
                           sec.axis = sec_axis(trans = ~. * scale_trans, name=element_blank(),
                                               labels = scales::number_format(accuracy=0.01)))
    
    # Cumulative Farm-CPU distribution
    if(nrow(farm.all) > 0){
        cum_color='gray35'
        for(mychrom in rownames(offsets)){
            
            # Get cumulative RMIP values
            mypoints = data.frame(Chr=mychrom, Pos=seq(0, offsets[mychrom, "length"], by=args$step), cum_rmip=0)
            for(i in 1:nrow(mypoints)){
                mypos = mypoints[i, "Pos"]
                nearby_hits = subset(farm.all, farm.all$Chr==mychrom & 
                                         abs(farm.all$Pos - mypos) <= args$winsize/2)
                mypoints$cum_rmip[i] = sum(nearby_hits$rmip)
            }
            
            # Get translated Y values
            mypoints$rmip.plot = rescale(mypoints$cum_rmip, from=c(0, max(farm.all$rmip)), to=c(0, max(glm$log.p)))
            
            # Add to plot
            myplot = myplot +
                geom_line(data=mypoints, mapping=aes(y=rmip.plot), color=cum_color, size=0.35) 
        }
        
        # Add horizontal reference line
        ref_line_y = rescale(args$rmip_cutoff, from=c(0, max(farm.all$rmip)), to=c(0, max(glm$log.p)))
        myplot = myplot +
            geom_hline(yintercept=ref_line_y, color=cum_color, size=0.35, linetype="dashed")
    }
    
    # Significant FarmCPU hits by RMIP
    if(nrow(farm) > 0){
        farm$rmip.plot = rescale(farm$rmip, from=c(0, max(farm$rmip)), to=c(0, max(glm$log.p)))
        rmip_color="black"
        myplot = myplot +
            #geom_segment(data=farm, mapping=aes(xend=get_offset(Chr, Pos), y=0, yend=rmip.plot), color=rmip_color, size=0.5) +
            geom_point(data=farm, mapping=aes(y=rmip.plot), size=2, color=rmip_color, fill=rmip_color, shape=25) 
    }
    
    return(myplot)
}

