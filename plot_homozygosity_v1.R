# Load functions
source('quantsmooth_modified_functions.R')

# Prepare function for loading data and plotting onto ideograms
get_plotting_data <- function(file) {
	if (is.na(file)) stop(c("Problem with file",file))
	a <- read.table(file,stringsAsFactors=F)[,2:3] 	# HARD-CODED, consider allowing arguments for chr/pos
	# Data now consist of two column data frame, containing columns with chr and position
	a <- a[a[,1] != "X",]
	a <- a[a[,1] != "Y",]
	CHR <- a[,1]
	MapInfo <- a[,2]
	x <- convert_df(data.frame(CHR,MapInfo))
	return(x)
}

#####################################
args   <- commandArgs(trailingOnly=T)
famid  <- args[[1]]
dadid  <- args[[2]]
momid  <- args[[3]]
suffix   <- "filtered.txt"	# Hard-coded, may need to adjust depending on file name
dad      <- get_plotting_data(paste(dadid, suffix, sep="."))
mom      <- get_plotting_data(paste(momid, suffix, sep="."))

# Plot ideograms
plot_outpathfile <- paste("indiv_fams/", args[[1]],".png",sep="")	# output directory
png(file=plot_outpathfile,width=700,height=700,res=75)
CHR<-seq(1,22)
MapInfo<-lengthChromosome(CHR,"bases") # position on chromosome
chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE, organism="hsa", main=paste("Family",args[[1]],sep=' '))

# Plot parents
spacing<-0.1
colors <- c("blue","red")
points(dad[,2],dad[,1]+1*(spacing),pch=15,col=colors[1])
points(mom[,2],mom[,1]+2*(spacing),pch=15,col=colors[2])

# Plot children if they exist
num_c <- length(args)-3  # famid,dad,mom are always defined
for (i in 1:num_c) { 
	colors <- c(colors,grey(i/(num_c+1)))
	data_file <- paste(args[[i+3]], suffix, sep=".") 
	tmp <- get_plotting_data(data_file)
	points(tmp[,2],tmp[,1]+(i+2)*(spacing),pch=15,col=colors[i+2])
}

colornames<- c("Dad", "Mom", "Child1","Child2","Child3","Child4","Child5")
colornames<- colornames[1:length(colors)]
legend(x="bottom", fill=rev(colors), legend=colornames)

dev.off()
