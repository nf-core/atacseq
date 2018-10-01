#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--frip_files"), type="character", default=NULL, help="Comma-separated list of FRiP score files. Each file should contain one line with FRiP score. ", metavar="frip_files"),
										make_option(c("-s", "--sample_ids"), type="character", default=NULL, help="Comma-separated list of sample ids associated with FRiP files. Must be unique and in same order as FRiP files input.", metavar="sampleids"),
									  make_option(c("-o", "--outdir"), type="character", default='./', help="Output directory", metavar="path"),
									  make_option(c("-p", "--outprefix"), type="character", default='plot_frip', help="Output prefix", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$frip_files)){
		print_help(opt_parser)
		stop("At least one FRiP score file must be supplied", call.=FALSE)
}
if (is.null(opt$sample_ids)){
		print_help(opt_parser)
		stop("Please provide sample ids associated with FRiP files.", call.=FALSE)
}

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir,recursive=TRUE)
}

FripFiles <- unlist(strsplit(opt$frip_files,","))
SampleIDs <- unlist(strsplit(opt$sample_ids,","))
if (length(FripFiles) != length(SampleIDs)) {
		print_help(opt_parser)
		stop("Number of sample ids must equal number of FRiP score files.", call.=FALSE)
}

################################################
################################################
## READ IN DATA                               ##
################################################
################################################

plot.dat <- data.frame()
for (idx in 1:length(FripFiles)) {
		sampleid = SampleIDs[idx]
		score.dat <- read.table(FripFiles[idx], header=FALSE)
		score.dat$name <- sampleid
		colnames(score.dat) <- c('frip','name')
		plot.dat <- rbind(plot.dat,score.dat)
}

################################################
################################################
## PLOTS                											##
################################################
################################################

PlotFile <- file.path(opt$outdir,paste(opt$outprefix,".plot.pdf",sep=""))
pdf(PlotFile,height=6,width=3*length(unique(plot.dat$name)))

## FRIP SCORE PLOT
plot  <- ggplot(plot.dat, aes(x=name, y=frip)) +
				 geom_bar(stat="identity",aes(colour=name,fill=name), position = "dodge", width = 0.8, alpha = 0.3) +
				 xlab("") +
				 ylab("FRiP score") +
				 ylim(0,1) +
				 theme(legend.position="none",
							 panel.grid.major = element_blank(),
							 panel.grid.minor = element_blank(),
							 panel.background = element_blank(),
							 axis.text.y = element_text(colour="black"),
							 axis.text.x= element_text(colour="black",face="bold"),
							 axis.line.x = element_line(size = 1, colour = "black", linetype = "solid"),
							 axis.line.y = element_line(size = 1, colour = "black", linetype = "solid")) +
				geom_text(aes(label = frip, x = name, y = frip), position = position_dodge(width = 0.8), vjust = -0.6)
print(plot)

dev.off()

################################################
################################################
################################################
################################################
