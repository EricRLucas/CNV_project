library(stringr)
library(Biostrings)
library(HardyWeinberg)

load('CNV_analysis_phase2_for_Gste2.Rdata')

# Get the coordinates of all the genes in the cluster
gste1.coords <- c(28598871, 28598871+945)
gste2.coords <- c(28597652, 28597652+988)
gste3.coords <- c(28601368, 28601368+912)
gste4.coords <- c(28595948, 28595948+920)
gste5.coords <- c(28594993, 28594993+809)
gste6.coords <- c(28593714, 28593714+1048)
gste7.coords <- c(28600501, 28600501+819)
gstu4.coords <- c(28591663, 28591663+812)

gst.region.coord <- c(28570000, 28620000) 
gst.region.plotting.indices <- c(min(which(counts.list[[1]]$Position >= gst.region.coord[1])), max(which(counts.list[[1]]$Position <= gst.region.coord[2])))

# Identify the samples with a detected duplication in the various genes in the region
encompass.duplication.detected.in.gste1 <- encompass.output.by.gene[['AGAP009195']]
encompass.duplication.detected.in.gste2 <- encompass.output.by.gene[['AGAP009194']]
encompass.duplication.detected.in.gste3 <- encompass.output.by.gene[['AGAP009197']]
encompass.duplication.detected.in.gste4 <- encompass.output.by.gene[['AGAP009193']]
encompass.duplication.detected.in.gste5 <- encompass.output.by.gene[['AGAP009192']]
encompass.duplication.detected.in.gste6 <- encompass.output.by.gene[['AGAP009191']]
encompass.duplication.detected.in.gste7 <- encompass.output.by.gene[['AGAP009196']]
encompass.duplication.detected.in.gstu4 <- encompass.output.by.gene[['AGAP009190']]
encompass.duplication.detected.in.gst.region <- unique(c(encompass.duplication.detected.in.gste1, 
                                                         encompass.duplication.detected.in.gste2, 
                                                         encompass.duplication.detected.in.gste3, 
                                                         encompass.duplication.detected.in.gste4, 
                                                         encompass.duplication.detected.in.gste5, 
                                                         encompass.duplication.detected.in.gste6, 
                                                         encompass.duplication.detected.in.gste7,  
                                                         encompass.duplication.detected.in.gstu4))

# Load table containing geographical id
meta <- read.table('/home/eric/Liverpool/phase2.AR1/samples/samples.meta.txt', header = T, row.names = 1, sep = '\t', quote = "", comment.char = "")
rownames(meta) <- sub('-', '_', rownames(meta))
# Create booleans vector that tell us whether each individual in the samples table has a duplication in one
# of the genes in the region
hasencomp <- rownames(meta) %in% encompass.duplication.detected.in.gst.region
# add these data to the meta table
meta$hasencomp = hasencomp

# Now, within each population, we want to know the proportion of individuals carrying a duplication
propencomp <- tapply(meta$hasencomp, meta$population, mean)
pop.sizes <- tapply(meta$population, meta$population, length)
propdup.for.plot <- rbind(propencomp, 1 - propencomp)

x11(width = 15)
par(mfrow = c(1,1), mar = c(2,4,4,8), oma = c(3,0,0,0), xpd = F)
barplot(propdup.for.plot, legend.text = c('Duplication', 'No duplication'), ylab = 'Proportion', names.arg = rep('', ncol(propdup.for.plot)), main = 'GSTE', args.legend = list(x=22.2, y=0.9))
mtext('Angola', side = 1, line = -0.1, at = 0.7, padj=1, cex=0.9)
mtext('M', side = 1, line = 1.9, at = 0.7, padj=1)
mtext('Burkina\nFaso', side = 1, line = -0.1, at = 1.9, padj=1, cex=0.9)
mtext('M', side = 1, line = 1.9, at = 1.9, padj=1)
mtext('Burkina\nFaso', side = 1, line = -0.1, at = 3.1, padj=1, cex=0.9)
mtext('S', side = 1, line = 1.9, at = 3.1, padj=1)
mtext('Côte\nd\'Ivoire', side = 1, line = -0.1, at = 4.3, padj=1, cex=0.9)
mtext('M', side = 1, line = 1.9, at = 4.3, padj=1)
mtext('Camer-\noon', side = 1, line = -0.1, at = 5.5, padj=1, cex=0.9)
mtext('S', side = 1, line = 1.9, at = 5.5, padj=1)
mtext('Mayotte', side = 1, line = -0.1, at = 6.7, padj=1, cex=0.9)
mtext('S', side = 1, line = 1.9, at = 6.7, padj=1)
mtext('Gabon', side = 1, line = -0.1, at = 7.9, padj=1, cex=0.9)
mtext('S', side = 1, line = 1.9, at = 7.9, padj=1)
mtext('Ghana', side = 1, line = -0.1, at = 9.1, padj=1, cex=0.9)
mtext('M', side = 1, line = 1.9, at = 9.1, padj=1)
mtext('Ghana', side = 1, line = -0.1, at = 10.3, padj=1, cex=0.9)
mtext('S', side = 1, line = 1.9, at = 10.3, padj=1)
mtext('The\nGambia', side = 1, line = -0.1, at = 11.5, padj=1, cex=0.9)
mtext('?', side = 1, line = 1.9, at = 11.5, padj=1)
mtext('Guinea', side = 1, line = -0.1, at = 12.7, padj=1, cex=0.9)
mtext('M', side = 1, line = 1.9, at = 12.7, padj=1)
mtext('Guinea', side = 1, line = -0.1, at = 13.9, padj=1, cex=0.9)
mtext('S', side = 1, line = 1.9, at = 13.9, padj=1)
mtext('Equ.\nGuinea', side = 1, line = -0.1, at = 15.1, padj=1, cex=0.9)
mtext('S', side = 1, line = 1.9, at = 15.1, padj=1)
mtext('Guinea\nBissau', side = 1, line = -0.1, at = 16.3, padj=1, cex=0.9)
mtext('M~S', side = 1, line = 1.9, at = 16.3, padj=1)
mtext('Kenya', side = 1, line = -0.1, at = 17.5, padj=1, cex=0.9)
mtext('S', side = 1, line = 1.9, at = 17.5, padj=1)
mtext('Uganda', side = 1, line = -0.1, at = 18.7, padj=1, cex=0.9)
mtext('S', side = 1, line = 1.9, at = 18.7, padj=1)
for (i in 1:length(pop.sizes)){
	mtext(paste('(', pop.sizes[i], ')', sep = ''), side = 1, line = 3, at = (-0.5) + 1.2*i, padj = 1)
}

# Now load up the results of the discordant read analysis
cat('Loading SSFA data\n')
SSFA.folder <- paste('/home/eric/Liverpool/CNV_/SSFA/v3_', chrom, '/phase2/GST_region', sep = '')
SSFA.all.files <- list.files(SSFA.folder)
SSFA.files <- SSFA.all.files[grepl('_output.csv$', SSFA.all.files)]
# Load each table 
FA.list <- list()
SS.list <- list()
FM.list <- list()
crosschrom.list <- list()
for (this.file in SSFA.files){
	# Extract the sample name information
	sample.name <- regmatches(this.file, regexpr('^A.\\d{4}_C.?(?=_)', this.file, perl=T))
	# Load the file and add the table to the counts.list
   	pos.table <- read.table(paste(SSFA.folder, this.file, sep='/'), header = T, sep='\t', stringsAsFactors = F)
	allpos.FA <- pos.table[grepl('FA.*>= 10', pos.table$Type), c('Position', 'Mate.position')]
	allpos.FA$Mate.position <- as.numeric(allpos.FA$Mate.position)
	allpos.FA <- allpos.FA[order(allpos.FA$Position),]
	allpos.SS <- pos.table[grepl('SS.*>= 10', pos.table$Type), c('Position', 'Mate.position')]
	allpos.SS$Mate.position <- as.numeric(allpos.SS$Mate.position)
	allpos.SS <- allpos.SS[order(allpos.SS$Position),]
	allpos.FM <- pos.table[grepl('FM.*>= 10', pos.table$Type), c('Position', 'Mate.position')]
	allpos.FM$Mate.position <- as.numeric(allpos.FM$Mate.position)
	allpos.FM <- allpos.FM[order(allpos.FM$Position),]
	allpos.crosschrom <- pos.table[grepl('crosschrom.*>= 10', pos.table$Type), c('Position', 'Mate.position')]
	allpos.crosschrom <- allpos.crosschrom[order(allpos.crosschrom$Position),]
	FA.list[[sample.name]] <- allpos.FA
	SS.list[[sample.name]] <- allpos.SS
	FM.list[[sample.name]] <- allpos.FM
	crosschrom.list[[sample.name]] <- allpos.crosschrom
}

# Get the faceaway reads in the gst region
FA.gst.list <- list()
for (i in 1:length(FA.list)){
	this.list <- FA.list[[i]]
	this.name <- names(FA.list)[i]
	# Use the following line if you want reads that start OR end in the gst region
	which.in.gst.FA <- ((this.list$Position >= gst.region.coord[1]) & (this.list$Position <= gst.region.coord[2])) | ((this.list$Mate.position >= gst.region.coord[1]) & (this.list$Mate.position <= gst.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the gst region
#	which.in.gst.FA <- ((this.list$Position >= gst.region.coord[1]) & (this.list$Mate.position <= gst.region.coord[2]))
	FA.gst.list[[this.name]] <- this.list[which.in.gst.FA, ]
}
SS.gst.list <- list()
for (i in 1:length(SS.list)){
	this.list <- SS.list[[i]]
	this.name <- names(SS.list)[i]
	# Use the following line if you want reads that start OR end in the gst region
	which.in.gst.SS <- ((this.list$Position >= gst.region.coord[1]) & (this.list$Position <= gst.region.coord[2])) | ((this.list$Mate.position >= gst.region.coord[1]) & (this.list$Mate.position <= gst.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the gst region
#	which.in.gst.SS <- ((this.list$Position >= gst.region.coord[1]) & (this.list$Mate.position <= gst.region.coord[2]))
	SS.gst.list[[this.name]] <- this.list[which.in.gst.SS, ]
}
FM.gst.list <- list()
for (i in 1:length(FM.list)){
	this.list <- FM.list[[i]]
	this.name <- names(FM.list)[i]
	# Use the following line if you want reads that start OR end in the gst region
	which.in.gst.FM <- ((this.list$Position >= gst.region.coord[1]) & (this.list$Position <= gst.region.coord[2])) | ((this.list$Mate.position >= gst.region.coord[1]) & (this.list$Mate.position <= gst.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the gst region
#	which.in.gst.FM <- ((this.list$Position >= gst.region.coord[1]) & (this.list$Mate.position <= gst.region.coord[2]))
	FM.gst.list[[this.name]] <- this.list[which.in.gst.FM, ]
}
crosschrom.gst.list <- list()
for (i in 1:length(crosschrom.list)){
	this.list <- crosschrom.list[[i]]
	this.name <- names(crosschrom.list)[i]
	which.in.gst.crosschrom <- (this.list$Position >= gst.region.coord[1]) & (this.list$Position <= gst.region.coord[2])
	crosschrom.gst.list[[this.name]] <- this.list[which.in.gst.crosschrom, ]
}

# Load up the results of the breakpoint detection 
cat('Loading GST region breakpoints data\n')
bp.folder <- paste('/home/eric/Liverpool/CNV_/breakpoints/' , chrom, '/phase2/GST_region', sep = '')
bp.all.files <- list.files(bp.folder)
bp.files <- bp.all.files[grepl('.csv', bp.all.files)]
# Load each table 
clipping.start.point.list <- list()
clipping.end.point.list <- list()
for (this.file in bp.files){
	# Extract the sample name information
	sample.name <- regmatches(this.file, regexpr('^A.\\d{4}_C.?(?=_)', this.file, perl=T))
	# Load the file and add the table to the counts.list
   	pos.table <- read.table(paste(bp.folder, this.file, sep='/'), header = T, sep='\t', stringsAsFactors = F)
	allpos.clipping.start.point <- pos.table[grepl('soft clipping start point mapq >= 10', pos.table$Type), c('Position', 'Clipped_sequence')]
	allpos.clipping.end.point <- pos.table[grepl('soft clipping end point mapq >= 10', pos.table$Type), c('Position', 'Clipped_sequence')]
	allpos.clipping.start.point <- allpos.clipping.start.point[order(allpos.clipping.start.point$Position), ]
	allpos.clipping.end.point <- allpos.clipping.end.point[order(allpos.clipping.end.point$Position), ]
	clipping.start.point.list[[sample.name]] <- allpos.clipping.start.point
	clipping.end.point.list[[sample.name]] <- allpos.clipping.end.point
}
clipping.start.point.gst.list <- list()
for (i in 1:length(clipping.start.point.list)){
	this.list <- clipping.start.point.list[[i]]
	this.name <- names(clipping.start.point.list)[i]
	which.in.gst.clipping.start.point <- (this.list$Position >= gst.region.coord[1]) & (this.list$Position <= gst.region.coord[2])
	this.initial.list <- this.list[which.in.gst.clipping.start.point, , drop = F]
	# The following for lines allow breakpoints to be kept only if their frequency exceeds a certain threshold. They 
	# are effectively pointless unless you set the threshold below above 0.
	reads.per.start.point <- tapply(this.initial.list$Position, as.factor(this.initial.list$Position), length)
	positions.to.keep <- this.initial.list$Position %in% names(reads.per.start.point)[reads.per.start.point > 0]
	clipping.start.point.gst.list[[this.name]] <- this.initial.list[positions.to.keep, ]
}
clipping.end.point.gst.list <- list()
for (i in 1:length(clipping.end.point.list)){
	this.list <- clipping.end.point.list[[i]]
	this.name <- names(clipping.end.point.list)[i]
	which.in.gst.clipping.end.point <- (this.list$Position >= gst.region.coord[1]) & (this.list$Position <= gst.region.coord[2])
	this.initial.list <- this.list[which.in.gst.clipping.end.point, , drop = F]
	reads.per.end.point <- tapply(this.initial.list$Position, as.factor(this.initial.list$Position), length)
	positions.to.keep <- this.initial.list$Position %in% names(reads.per.end.point)[reads.per.end.point > 0]
	clipping.end.point.gst.list[[this.name]] <- this.initial.list[positions.to.keep, ]
}
cat('Loading 5595 region breakpoints data\n')
bp.5595.folder <- paste('/home/eric/Liverpool/CNV_/breakpoints/', chrom, '/phase2/5595_region', sep = '')
bp.5595.all.files <- list.files(bp.5595.folder)
bp.5595.files <- bp.5595.all.files[grepl('.csv', bp.5595.all.files)]
# Load each table 
clipping.start.point.5595.list <- list()
clipping.end.point.5595.list <- list()
for (this.file in bp.5595.files){
	# Extract the sample name information
	sample.name <- regmatches(this.file, regexpr('^A.\\d{4}_C.?(?=_)', this.file, perl=T))
	# Load the file and add the table to the counts.list
   	pos.table <- read.table(paste(bp.5595.folder, this.file, sep='/'), header = T, sep='\t', stringsAsFactors = F)
	allpos.clipping.start.point <- pos.table[grepl('soft clipping start point mapq >= 10', pos.table$Type), c('Position', 'Clipped_sequence')]
	allpos.clipping.end.point <- pos.table[grepl('soft clipping end point mapq >= 10', pos.table$Type), c('Position', 'Clipped_sequence')]
	allpos.clipping.start.point <- allpos.clipping.start.point[order(allpos.clipping.start.point$Position), ]
	allpos.clipping.end.point <- allpos.clipping.end.point[order(allpos.clipping.end.point$Position), ]
	clipping.start.point.5595.list[[sample.name]] <- allpos.clipping.start.point
	clipping.end.point.5595.list[[sample.name]] <- allpos.clipping.end.point
}

# Create a function that can be used to plot discordant read pairs
add.faceaways <- function(pair.matrix, this.col = 'red', yrange = c(0,1)){
	jitter.values <- seq(yrange[1], yrange[2], length.out = nrow(pair.matrix))
	y.values <- cbind(jitter.values, jitter.values)
	points(as.matrix(pair.matrix), y.values, col = this.col)
	for (i in 1:nrow(pair.matrix)){
		lines(pair.matrix[i,], y.values[i,], col = this.col)
	}
}

# Create a function that can be used to plot crosschrom and breakpoint values
add.breakpoints <- function(points.vector, this.col = 'red', yrange = c(0,1)){
	jitter.values <- seq(yrange[1], yrange[2], length.out = length(points.vector))
	points(points.vector, jitter.values, col = this.col, pch = 10)
}

# Get the mapq0 proportions:
mapq0.proportions.allchrom <- read.table('/home/eric/Liverpool/CNV_v2/mapq_proportions_allchrom.txt', header = T, sep = '\t', check.names = F)
mapq0.proportions <- subset(mapq0.proportions.allchrom, Chrom == chrom)[,2:4]
# Restrict the table to the region studied here 
mapq0.proportions <- subset(mapq0.proportions, (Position >= gst.region.coord[1]) & (Position <= gst.region.coord[2]))
total.coverage <- mapq0.proportions$'Count mapq > 0' + mapq0.proportions$'Count mapq = 0'
mapq0.proportions$mapq0.prop <- mapq0.proportions$'Count mapq = 0' / total.coverage
# In some cases, there is no coverage, so we get infinite values. Let's deal with those
mapq0.proportions$mapq0.prop[total.coverage == 0] <- 0

# Create a list of colours that will be used for plotting the duplication types
duplication.colours <- c(rgb(0.6509804, 0.8078431, 0.8901961), rgb(0.1215686, 0.4705882, 0.7058824), rgb(0.6980392, 0.8745098, 0.5411765), rgb(0.2, 0.627451, 0.172549), rgb(0.9843137, 0.6039216, 0.6), rgb(0.8901961, 0.1019608, 0.1098039), rgb(0.9921569, 0.7490196, 0.4352941), rgb(1, 0.4980392, 0), rgb(0.7921569, 0.6980392, 0.8392157), rgb(0.4156863, 0.2392157, 0.6039216), rgb(1,0,1,0.7))#, rgb(0.6941176, 0.3490196, 0.1568627), rgb(0.5,0.5,0,0.8), rgb(0,0.5,0.5,0.8), rgb(0.5, 0, 0.5, 0.8))

plot.gst <- function(this.sample, list.of.dups = NULL, FA = T, SS=T, FM=T, XC=T, BP=T, start.index = gst.region.plotting.indices[1], end.index = gst.region.plotting.indices[2]){
	this.sample.name <- regmatches(this.sample, regexpr('A.\\d\\d\\d\\d_C', this.sample))
	plot(counts.list[[this.sample]]$Position[start.index : end.index], counts.list[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		if (list.of.dups['Dup11'] >= 1){
			rect(Dup11.bp[1], -0.6, Dup11.bp[2], 0, col = rgb(1, 1, 0, 0.8), border = rgb(1, 1, 0, 0.8))
			text(Dup11.bp[2], -0.4, 'Dup11', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup10'] >= 1){
			rect(Dup10.bp[1], -0.6, Dup10.bp[2], 0, col = rgb(0, 1, 1, 0.8), border = rgb(0, 1, 1, 0.8))
			text(Dup10.bp[2], -0.4, 'Dup10', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup9'] >= 1){
			rect(Dup9.bp[1], -0.6, Dup9.bp[2], 0, col = rgb(0.5, 0, 0, 0.8), border = rgb(0.5, 0, 0, 0.8))
			text(Dup9.bp[2], -0.4, 'Dup9', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup8'] >= 1){
			rect(Dup8.bp[1], -0.6, Dup8.bp[2], 0, col = rgb(0, 0, 1, 0.8), border = rgb(0, 0, 1, 0.8))
			text(Dup8.bp[2], -0.4, 'Dup8', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup7'] >= 1){
			rect(Dup7.bp[1], -0.6, Dup7.bp[2], 0, col = rgb(1, 0, 0, 0.8), border = rgb(1, 0, 0, 0.8))
			text(Dup7.bp[2], -0.4, 'Dup7', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup6'] >= 1){
			rect(Dup6.bp[1], -0.6, Dup6.bp[2], 0, col = rgb(0.5, 0.5, 0, 0.8), border = rgb(0.5, 0.5, 0, 0.8))
			text(Dup6.bp[2], -0.4, 'Dup6', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup5'] >= 1){
			rect(Dup5.bp[1], -0.6, Dup5.bp[2], 0, col = rgb(0.5, 0, 0.5, 0.8), border = rgb(0.5, 0, 0.5, 0.8))
			text(Dup5.bp[2], -0.4, 'Dup5', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup4'] >= 1){
			rect(Dup4.bp[1], -0.6, Dup4.bp[2], 0, col = rgb(0, 0, 0.5, 0.8), border = rgb(0, 0, 0.5, 0.8))
			text(Dup4.bp[2], -0.4, 'Dup4', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup3'] >= 1){
			rect(Dup3.bp[1], -0.6, Dup3.bp[2], 0, col = rgb(1, 0, 1, 0.8), border = rgb(1, 0, 1, 0.8))
			text(Dup3.bp[2], -0.4, 'Dup3', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup2'] >= 1){
			rect(Dup2.bp[1], -0.6, Dup2.bp[2], 0, col = rgb(0, 0.5, 0, 0.8), border = rgb(0, 0.5, 0, 0.8))
			text(Dup2.bp[2], -0.4, 'Dup2', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup1'] >= 1){
			rect(Dup1.bp[1], -0.6, Dup1.bp[2], 0, col = rgb(0, 1, 0, 0.8), border = rgb(0, 1, 0, 0.8))
			text(Dup1.bp[2], -0.4, 'Dup1', cex = 1.2, adj = c(1,0))
		}
	}
	if (FA){
		if (nrow(FA.gst.list[[this.sample.name]]) > 0)
			add.faceaways(FA.gst.list[[this.sample.name]], this.col = 'blue', yrange = c(0,3))
	}
	if (SS){
		if (nrow(SS.gst.list[[this.sample.name]]) > 0)
			add.faceaways(SS.gst.list[[this.sample.name]], this.col = 'cyan', yrange = c(3,6))
	}
	if (FM){
		if (nrow(FM.gst.list[[this.sample.name]]) > 0)
			add.faceaways(FM.gst.list[[this.sample.name]], this.col = 'green', yrange = c(6,9))
	}
	if (XC){
		if (nrow(crosschrom.gst.list[[this.sample.name]])) 
			add.breakpoints(crosschrom.gst.list[[this.sample.name]]$Position, this.col = 'red', yrange = c(9,12))
	}
	if (BP){
		if (nrow(clipping.end.point.gst.list[[this.sample.name]]))
			add.breakpoints(clipping.end.point.gst.list[[this.sample.name]][, 'Position'], this.col = 'brown', yrange = c(0,12))
		if (nrow(clipping.start.point.gst.list[[this.sample.name]]))
			add.breakpoints(clipping.start.point.gst.list[[this.sample.name]][, 'Position'], this.col = 'pink', yrange = c(0,12))
	}
	these.cnv.states <- counts.list[[this.sample]]$CNV[start.index : end.index]
	lines(counts.list[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2)
	abline(v = gste1.coords[1])
	abline(v = gste1.coords[2])
	text(mean(gste1.coords), 9, 'GSTE1', srt=90, adj = 0, col = 'black')
	abline(v = gste2.coords[1], col = 'purple')
	abline(v = gste2.coords[2], col = 'purple')
	text(mean(gste2.coords), 9, 'GSTE2', srt=90, adj = 0, col = 'purple')
	abline(v = gste3.coords[1], col = 'orange')
	abline(v = gste3.coords[2], col = 'orange')
	text(mean(gste3.coords), 9, 'GSTE3', srt=90, adj = 0, col = 'orange')
	abline(v = gste4.coords[1], col = 'magenta')
	abline(v = gste4.coords[2], col = 'magenta')
	text(mean(gste4.coords), 9, 'GSTE4', srt=90, adj = 0, col = 'magenta')
	abline(v = gste5.coords[1], col = 'brown')
	abline(v = gste5.coords[2], col = 'brown')
	text(mean(gste5.coords), 9, 'GSTE5', srt=90, adj = 0, col = 'brown')
	abline(v = gste6.coords[1], col = 'green')
	abline(v = gste6.coords[2], col = 'green')
	text(mean(gste6.coords), 9, 'GSTE6', srt=90, adj = 0, col = 'green')
	abline(v = gste7.coords[1], col = 'violet')
	abline(v = gste7.coords[2], col = 'violet')
	text(mean(gste7.coords), 9, 'GSTE7', srt=90, adj = 0, col = 'violet')
	abline(v = gstu4.coords[1], col = 'red')
	abline(v = gstu4.coords[2], col = 'red')
	text(mean(gstu4.coords), 9, 'GSTU4', srt=90, adj = 0, col = 'red')
}


plot.all.gst <- function(list.of.all, list.of.dupl = encompass.duplication.detected.in.gst.region, FA=T, SS=T, FM=T, XC=T, BP=T, start.index = gst.region.plotting.indices[1], end.index = gst.region.plotting.indices[2]){
	x11()
	list.of.all <- sub('-', '_', list.of.all)
	list.of.dupl <- sub('-', '_', list.of.dupl)
	for (this.sample in list.of.all){
		if (this.sample %in% list.of.dupl)
			par(bg = 'beige')
		else
			par(bg = 'white')
		plot.gst(this.sample, NULL, FA, SS, FM, XC, BP, start.index, end.index)
		locator(1)
	}
}

plot.all.gst.with.read.dups <- function(list.of.samples, matrix.of.read.dups = read.based.gst.duplications, FA = T, SS = T, FM = T, XC=T, BP=T, start.index = gst.region.plotting.indices[1], end.index = gst.region.plotting.indices[2]){
	x11()
	x.midpoint <- mean(c(counts.list[[1]]$Position[start.index], counts.list[[1]]$Position[end.index]))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- sub('-', '_', list.of.samples[i])
		plot.gst(this.sample, matrix.of.read.dups[this.sample,], FA, SS, FM, XC, BP, start.index, end.index)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

# Use the following definitions for the different duplications (numbers represent start range (1st row) and end range
# (2nd row) of the read pairs). These values were determined after visual inspection of the data. 
Dup1.fa <- data.frame(c(28596750, 28598550), c(28597050, 28598850))
Dup1.bp <- c(28596818, 28598850)
Dup1.bp.seq <- c('TTTTG', 'CGTTT')
#
# The following definition includes some false positives. 
Dup2.bp <- c(28596390, 28598923)
Dup2.bp.seq <- c('GGGGG', 'TTCCC')
#
Dup3.fa <- data.frame(c(28590500, 28592950), c(28590800, 28593250))
Dup3.bp <- c(28590597, 28593254)
Dup3.bp.seq <- c('TCAAA', 'AGGGC')
#
Dup4.fa <- data.frame(c(28595050, 28598750), c(28595350, 28599050))
Dup4.bp <- c(28595162, 28599081)
Dup4.bp.seq <- c('TTCTA', 'AGAAC')
#
Dup5.bp <- c(28593122, 28598971)
Dup5.bp.seq <- c('GTCAT', 'ATTTA')
#
Dup6.fa <- data.frame(c(28596250, 28601900), c(28596550, 28602200))
Dup6.bp <- c(28596241, 28602177)
Dup6.bp.seq <- c('ACAAC', 'GAAGC')
#
Dup7.xc.start <- data.frame(c(28597400, 3696450), c(28597700, 3696750), c('3R', '2L'))
Dup7.xc.end <- data.frame(c(28603950, 26597300), c(28604250, 26597600), c('3R', 'UNKN'))
Dup7.bp <- c(28597504, 28604250)
Dup7.bp.seq <- c('GTCCA', 'GCTGT')
#
Dup8.bp <- c(28594797, 28602349)
Dup8.bp.seq <- c('GTCCC', 'CAGGG')
#
Dup9.fa <- data.frame(c(28591050, 28600850), c(28591350, 28601150))
Dup9.bp <- c(28591140, 28601188)
Dup9.bp.seq <- c('AGAAG', 'GATGA')
#
Dup10.fa <- data.frame(c(28593550, 28603350), c(28593850, 28603650))
Dup10.bp <- c(28593642, 28603786)
Dup10.bp.seq <- c('TCGCT', 'AAGAC')
#
Dup11.xc.start <- data.frame(c(28581250, 29210650), c(28581550, 29210950), c('3R', 'UNKN'))
Dup11.xc.end <- data.frame(c(28604650, 29210650), c(28604950, 29210950), c('3R', 'UNKN'))
Dup11.bp <- c(28581256, 28604994)
Dup11.bp.seq <- c('CCATT', 'GGTAA')

# We will store the duplication type information in a matrix
read.based.gst.duplications <- matrix(0, length(counts.list), 12, dimnames = list(names(counts.list), c('Dup0', 'Dup1', 'Dup2', 'Dup3', 'Dup4', 'Dup5', 'Dup6', 'Dup7', 'Dup8', 'Dup9', 'Dup10', 'Dup11')))
# Now go through each individual and determine which duplications they have based on their diagnostic reads. 
for (this.sample in names(counts.list)){
	# get the FA and SS reads encompassed in the general region
	these.FA <- FA.gst.list[[this.sample]]
	these.XC <- crosschrom.gst.list[[this.sample]]
	these.XC$Matechrom <- sub(':.*', '', these.XC$Mate.position)
	these.XC$Matepos <- as.integer(sub('.*:' , '', these.XC$Mate.position))
	these.CSP <- clipping.start.point.gst.list[[this.sample]]
	these.CEP <- clipping.end.point.gst.list[[this.sample]]
	# Count the number of supporting reads
	num.Dup1.fa.reads <- sum((these.FA$Position > Dup1.fa[1,1]) & these.FA$Position < Dup1.fa[1,2]
	                       & (these.FA$Mate.position > Dup1.fa[2,1]) & (these.FA$Mate.position < Dup1.fa[2,2]))
	Dup1.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup1.bp[1]]
	num.Dup1.start.reads <- sum(substr(reverse(Dup1.CEP.seq), 1, 5) == Dup1.bp.seq[1])
	Dup1.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup1.bp[2]]
	num.Dup1.end.reads <- sum(substr(Dup1.CSP.seq, 1, 5) == Dup1.bp.seq[2])
	num.Dup1.reads <- num.Dup1.fa.reads + num.Dup1.start.reads + num.Dup1.end.reads
	#
	Dup2.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup2.bp[1]]
	num.Dup2.start.reads <- sum(substr(reverse(Dup2.CEP.seq), 1, 5) == Dup2.bp.seq[1])
	Dup2.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup2.bp[2]]
	num.Dup2.end.reads <- sum(substr(Dup2.CSP.seq, 1, 5) == Dup2.bp.seq[2])
	num.Dup2.reads <- min(num.Dup2.start.reads, num.Dup2.end.reads)
	#
	num.Dup3.fa.reads <- sum((these.FA$Position > Dup3.fa[1,1]) & these.FA$Position < Dup3.fa[1,2]
	                       & (these.FA$Mate.position > Dup3.fa[2,1]) & (these.FA$Mate.position < Dup3.fa[2,2]))
	Dup3.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup3.bp[1]]
	num.Dup3.start.reads <- sum(substr(reverse(Dup3.CEP.seq), 1, 5) == Dup3.bp.seq[1])
	Dup3.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup3.bp[2]]
	num.Dup3.end.reads <- sum(substr(Dup3.CSP.seq, 1, 5) == Dup3.bp.seq[2])
	num.Dup3.reads <- num.Dup3.fa.reads + num.Dup3.start.reads + num.Dup3.end.reads
	#
	num.Dup4.fa.reads <- sum((these.FA$Position > Dup4.fa[1,1]) & these.FA$Position < Dup4.fa[1,2]
	                       & (these.FA$Mate.position > Dup4.fa[2,1]) & (these.FA$Mate.position < Dup4.fa[2,2]))
	Dup4.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup4.bp[1]]
	num.Dup4.start.reads <- sum(substr(reverse(Dup4.CEP.seq), 1, 5) == Dup4.bp.seq[1])
	Dup4.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup4.bp[2]]
	num.Dup4.end.reads <- sum(substr(Dup4.CSP.seq, 1, 5) == Dup4.bp.seq[2])
	num.Dup4.reads <- num.Dup4.fa.reads + num.Dup4.start.reads + num.Dup4.end.reads
	#
	Dup5.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup5.bp[1]]
	num.Dup5.start.reads <- sum(substr(reverse(Dup5.CEP.seq), 1, 5) == Dup5.bp.seq[1])
	Dup5.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup5.bp[2]]
	num.Dup5.end.reads <- sum(substr(Dup5.CSP.seq, 1, 5) == Dup5.bp.seq[2])
	num.Dup5.reads <- num.Dup5.start.reads + num.Dup5.end.reads
	#
	num.Dup6.fa.reads <- sum((these.FA$Position > Dup6.fa[1,1]) & these.FA$Position < Dup6.fa[1,2]
	                        & (these.FA$Mate.position > Dup6.fa[2,1]) & (these.FA$Mate.position < Dup6.fa[2,2]))
	Dup6.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup6.bp[1]]
	num.Dup6.start.reads <- sum(substr(reverse(Dup6.CEP.seq), 1, 5) == Dup6.bp.seq[1])
	Dup6.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup6.bp[2]]
	num.Dup6.end.reads <- sum(substr(Dup6.CSP.seq, 1, 5) == Dup6.bp.seq[2])
	num.Dup6.reads <- num.Dup6.fa.reads + num.Dup6.start.reads + num.Dup6.end.reads
	#
	num.Dup7.xc.start.reads <- sum((these.XC$Position > Dup7.xc.start[1,1]) & (these.XC$Position < Dup7.xc.start[1,2])
	                             & (these.XC$Matechrom == Dup7.xc.start[2,3]) & (these.XC$Matepos > Dup7.xc.start[2,1]) & (these.XC$Matepos < Dup7.xc.start[2,2]), na.rm = T)
	num.Dup7.xc.end.reads <- sum((these.XC$Position > Dup7.xc.end[1,1]) & (these.XC$Position < Dup7.xc.end[1,2])
	                             & (these.XC$Matechrom == Dup7.xc.end[2,3]) & (these.XC$Matepos > Dup7.xc.end[2,1]) & (these.XC$Matepos < Dup7.xc.end[2,2]), na.rm = T)
	Dup7.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup7.bp[1]]
	num.Dup7.start.reads <- sum(substr(reverse(Dup7.CEP.seq), 1, 5) == Dup7.bp.seq[1])
	Dup7.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup7.bp[2]]
	num.Dup7.end.reads <- sum(substr(Dup7.CSP.seq, 1, 5) == Dup7.bp.seq[2])
	num.Dup7.reads <- num.Dup7.xc.start.reads + num.Dup7.xc.end.reads + num.Dup7.start.reads + num.Dup7.end.reads
	# 
	Dup8.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup8.bp[1]]
	num.Dup8.start.reads <- sum(substr(reverse(Dup8.CEP.seq), 1, 5) == Dup8.bp.seq[1])
	Dup8.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup8.bp[2]]
	num.Dup8.end.reads <- sum(substr(Dup8.CSP.seq, 1, 5) == Dup8.bp.seq[2])
	num.Dup8.reads <- num.Dup8.start.reads + num.Dup8.end.reads
	#
	num.Dup9.fa.reads <- sum((these.FA$Position > Dup9.fa[1,1]) & these.FA$Position < Dup9.fa[1,2]
	                       & (these.FA$Mate.position > Dup9.fa[2,1]) & (these.FA$Mate.position < Dup9.fa[2,2]))
	Dup9.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup9.bp[1]]
	num.Dup9.start.reads <- sum(substr(reverse(Dup9.CEP.seq), 1, 5) == Dup9.bp.seq[1])
	Dup9.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup9.bp[2]]
	num.Dup9.end.reads <- sum(substr(Dup9.CSP.seq, 1, 5) == Dup9.bp.seq[2])
	num.Dup9.reads <- num.Dup9.fa.reads + num.Dup9.start.reads + num.Dup9.end.reads
	#
	num.Dup10.fa.reads <- sum((these.FA$Position > Dup10.fa[1,1]) & these.FA$Position < Dup10.fa[1,2]
	                       & (these.FA$Mate.position > Dup10.fa[2,1]) & (these.FA$Mate.position < Dup10.fa[2,2]))
	Dup10.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup10.bp[1]]
	num.Dup10.start.reads <- sum(substr(reverse(Dup10.CEP.seq), 1, 5) == Dup10.bp.seq[1])
	Dup10.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup10.bp[2]]
	num.Dup10.end.reads <- sum(substr(Dup10.CSP.seq, 1, 5) == Dup10.bp.seq[2])
	num.Dup10.reads <- num.Dup10.fa.reads + num.Dup10.start.reads + num.Dup10.end.reads
	#
	num.Dup11.xc.start.reads <- sum((these.XC$Position > Dup11.xc.start[1,1]) & (these.XC$Position < Dup11.xc.start[1,2])
	                             & (these.XC$Matechrom == Dup11.xc.start[2,3]) & (these.XC$Matepos > Dup11.xc.start[2,1]) & (these.XC$Matepos < Dup11.xc.start[2,2]), na.rm = T)
	num.Dup11.xc.end.reads <- sum((these.XC$Position > Dup11.xc.end[1,1]) & (these.XC$Position < Dup11.xc.end[1,2])
	                             & (these.XC$Matechrom == Dup11.xc.end[2,3]) & (these.XC$Matepos > Dup11.xc.end[2,1]) & (these.XC$Matepos < Dup11.xc.end[2,2]), na.rm = T)
	Dup11.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup11.bp[1]]
	num.Dup11.start.reads <- sum(substr(reverse(Dup11.CEP.seq), 1, 5) == Dup11.bp.seq[1])
	Dup11.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup11.bp[2]]
	num.Dup11.end.reads <- sum(substr(Dup11.CSP.seq, 1, 5) == Dup11.bp.seq[2])
	num.Dup11.reads <- num.Dup11.xc.start.reads + num.Dup11.xc.end.reads + num.Dup11.start.reads + num.Dup11.end.reads
	# For each of these, if there are at least 2 supporting reads, score it as a duplication. If there is only one
	# supporting read, score it tentatively
	read.based.gst.duplications[this.sample, 2:12] <- as.logical(c(num.Dup1.reads, num.Dup2.reads, num.Dup3.reads, num.Dup4.reads, num.Dup5.reads, num.Dup6.reads, num.Dup7.reads, num.Dup8.reads, num.Dup9.reads, num.Dup10.reads, num.Dup11.reads))
	if (num.Dup1.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup1'] <- 0.5
	if (num.Dup2.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup2'] <- 0.5
	if (num.Dup3.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup3'] <- 0.5
	if (num.Dup4.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup4'] <- 0.5
	if (num.Dup5.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup5'] <- 0.5
	if (num.Dup6.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup6'] <- 0.5
	if (num.Dup7.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup7'] <- 0.5
	if (num.Dup8.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup8'] <- 0.5
	if (num.Dup9.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup9'] <- 0.5
	if (num.Dup10.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup10'] <- 0.5
	if (num.Dup11.reads == 1)
		read.based.gst.duplications[this.sample, 'Dup11'] <- 0.5
	# If the sample has a duplication, but is none of the specific ones, we record it as Dup0
	if ((sum(read.based.gst.duplications[this.sample, ] >= 1) == 0) & (this.sample %in% encompass.duplication.detected.in.gst.region))
		read.based.gst.duplications[this.sample, 'Dup0'] <- 1
}

# Compare the duplications detected by coverage or reads
cov.based.duplications.names <- encompass.duplication.detected.in.gst.region
read.based.duplications.names <- rownames(read.based.gst.duplications)[as.logical(apply(read.based.gst.duplications[,2:12] >= 1, 1, sum))]

cov.based.negatives <- setdiff(read.based.duplications.names, cov.based.duplications.names)
read.based.negatives <- setdiff(cov.based.duplications.names, read.based.duplications.names)

# Output a table of each duplication type and the country they are found in
Dup0.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup0'] >= 1], 'population', drop = F]
Dup1.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup1'] >= 1], 'population', drop = F]
Dup2.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup2'] >= 1], 'population', drop = F]
Dup3.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup3'] >= 1], 'population', drop = F]
Dup4.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup4'] >= 1], 'population', drop = F]
Dup5.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup5'] >= 1], 'population', drop = F]
Dup6.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup6'] >= 1], 'population', drop = F]
Dup7.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup7'] >= 1], 'population', drop = F]
Dup8.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup8'] >= 1], 'population', drop = F]
Dup9.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup9'] >= 1], 'population', drop = F]
Dup10.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup10'] >= 1], 'population', drop = F]
Dup11.meta <- meta[rownames(read.based.gst.duplications)[read.based.gst.duplications[,'Dup11'] >= 1], 'population', drop = F]

# Let's output some lists of individuals carrying each duplication
write(rownames(Dup0.meta), file = 'Dup0_samples.txt', ncolumns = 1)
write(rownames(Dup1.meta), file = 'Dup1_samples.txt', ncolumns = 1)
write(rownames(Dup2.meta), file = 'Dup2_samples.txt', ncolumns = 1)
write(rownames(Dup3.meta), file = 'Dup3_samples.txt', ncolumns = 1)
write(rownames(Dup4.meta), file = 'Dup4_samples.txt', ncolumns = 1)
write(rownames(Dup5.meta), file = 'Dup5_samples.txt', ncolumns = 1)
write(rownames(Dup6.meta), file = 'Dup6_samples.txt', ncolumns = 1)
write(rownames(Dup7.meta), file = 'Dup7_samples.txt', ncolumns = 1)
write(rownames(Dup8.meta), file = 'Dup8_samples.txt', ncolumns = 1)
write(rownames(Dup9.meta), file = 'Dup9_samples.txt', ncolumns = 1)
write(rownames(Dup10.meta), file = 'Dup10_samples.txt', ncolumns = 1)
write(rownames(Dup11.meta), file = 'Dup11_samples.txt', ncolumns = 1)

Dup0.samples <- rownames(Dup0.meta)
Dup1.samples <- rownames(Dup1.meta)
Dup2.samples <- rownames(Dup2.meta)
Dup3.samples <- rownames(Dup3.meta)
Dup4.samples <- rownames(Dup4.meta)
Dup5.samples <- rownames(Dup5.meta)
Dup6.samples <- rownames(Dup6.meta)
Dup7.samples <- rownames(Dup7.meta)
Dup8.samples <- rownames(Dup8.meta)
Dup9.samples <- rownames(Dup9.meta)
Dup10.samples <- rownames(Dup10.meta)
Dup11.samples <- rownames(Dup11.meta)


######################################
# PRECISELY IDENTIFYING BREAK POINTS #
######################################

# Load the genome
genome <- readDNAStringSet('/home/eric/Liverpool/AR3/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa')

# Dup1:
assumed.Dup1.start.position <- 28596818
assumed.Dup1.end.position <- 28598850
all.Dup1.start.sequences <- character()
for (sample.name in sub('-', '_', Dup1.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup1.start.position, 'Clipped_sequence']
	all.Dup1.start.sequences <- c(all.Dup1.start.sequences, these.start.sequences)
}
longest.Dup1.start.sequence <- all.Dup1.start.sequences[which.max(nchar(all.Dup1.start.sequences))]
# Do the same at the end of the sequence
all.Dup1.end.sequences <- character()
for (sample.name in sub('-', '_', Dup1.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup1.end.position, 'Clipped_sequence']
	all.Dup1.end.sequences <- c(all.Dup1.end.sequences, these.end.sequences)
}
longest.Dup1.end.sequence <- all.Dup1.end.sequences[which.max(nchar(all.Dup1.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup1 <- substr(genome[4], assumed.Dup1.end.position-100, assumed.Dup1.end.position-1)
right.Dup1 <- substr(genome[4], assumed.Dup1.start.position+1, assumed.Dup1.start.position+100)
# The sequences are offset by one, this is because there is a T either side of the breakpoint, so we can't know
# for sure where the break happened, and the T aligns in both directions (ie: whichever side of the T is after the
# break, it still aligns because it's still a T
# We check that the clipped sequences are indeed as we expect. 
grepl(longest.Dup1.start.sequence, left.Dup1)
grepl(longest.Dup1.end.sequence, right.Dup1)
# The breakpoint looks like this (the T is present at both ends and therefore could be either side of the breakpoint)
# AGAAGCGAATTCCTGTTTT T CGTTTGAATGGCGTTTCGGGCT
#    end of the dup ^   ^ start of the dup
#                   ^   ^ 
#   position 28598848   position 28596820

# Dup2:
assumed.Dup2.start.position <- 28596390
assumed.Dup2.end.position <- 28598923
all.Dup2.start.sequences <- character()
for (sample.name in sub('-', '_', Dup2.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup2.start.position, 'Clipped_sequence']
	all.Dup2.start.sequences <- c(all.Dup2.start.sequences, these.start.sequences)
}
longest.Dup2.start.sequence <- all.Dup2.start.sequences[which.max(nchar(all.Dup2.start.sequences))]
# Do the same at the end of the sequence
all.Dup2.end.sequences <- character()
for (sample.name in sub('-', '_', Dup2.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup2.end.position, 'Clipped_sequence']
	all.Dup2.end.sequences <- c(all.Dup2.end.sequences, these.end.sequences)
}
longest.Dup2.end.sequence <- all.Dup2.end.sequences[which.max(nchar(all.Dup2.end.sequences))]
# Again, the clipped reads align in multiple places in the genome. 

# Dup3:
assumed.Dup3.start.position <- 28590597
assumed.Dup3.end.position <- 28593254
all.Dup3.start.sequences <- character()
for (sample.name in sub('-', '_', Dup3.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup3.start.position, 'Clipped_sequence']
	all.Dup3.start.sequences <- c(all.Dup3.start.sequences, these.start.sequences)
}
longest.Dup3.start.sequence <- all.Dup3.start.sequences[which.max(nchar(all.Dup3.start.sequences))]
# Do the same at the end of the sequence
all.Dup3.end.sequences <- character()
for (sample.name in sub('-', '_', Dup3.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup3.end.position, 'Clipped_sequence']
	all.Dup3.end.sequences <- c(all.Dup3.end.sequences, these.end.sequences)
}
longest.Dup3.end.sequence <- all.Dup3.end.sequences[which.max(nchar(all.Dup3.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point.
left.Dup3 <- substr(genome[4], assumed.Dup3.end.position-100, assumed.Dup3.end.position-1)
right.Dup3 <- substr(genome[4], assumed.Dup3.start.position+1, assumed.Dup3.start.position+100)
# The sequences are offset by one, this is because there is a T either side of the breakpoint, so we can't know
# for sure where the break happened, and the T aligns in both directions (ie: whichever side of the T is after the
# break, it still aligns because it's still a T
# We check that the clipped sequences are indeed as we expect. 
grepl(longest.Dup3.start.sequence, left.Dup3)
grepl(longest.Dup3.end.sequence, right.Dup3)
# Those two tests come out as false, but only because of a single SNP in each sequence. They are clearly matching
# correctly. 
# The breakpoint looks like this (the T is present at both ends and therefore could be either side of the breakpoint)
# GAAAAATCTAGTTGAAACT T AGGGCTATAGTATATATC
#    end of the dup ^   ^ start of the dup
#                   ^   ^ 
#   position 28593252   position 28590599

# Dup4:
assumed.Dup4.start.position <- 28595162
assumed.Dup4.end.position <- 28599081
all.Dup4.start.sequences <- character()
for (sample.name in sub('-', '_', Dup4.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup4.start.position, 'Clipped_sequence']
	all.Dup4.start.sequences <- c(all.Dup4.start.sequences, these.start.sequences)
}
longest.Dup4.start.sequence <- all.Dup4.start.sequences[which.max(nchar(all.Dup4.start.sequences))]
# Do the same at the end of the sequence
all.Dup4.end.sequences <- character()
for (sample.name in sub('-', '_', Dup4.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup4.end.position, 'Clipped_sequence']
	all.Dup4.end.sequences <- c(all.Dup4.end.sequences, these.end.sequences)
}
longest.Dup4.end.sequence <- all.Dup4.end.sequences[which.max(nchar(all.Dup4.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point.
left.Dup4 <- substr(genome[4], assumed.Dup4.end.position-100, assumed.Dup4.end.position-1)
right.Dup4 <- substr(genome[4], assumed.Dup4.start.position+1, assumed.Dup4.start.position+100)
# We check that the clipped sequences are indeed as we expect. 
grepl(longest.Dup4.start.sequence, left.Dup4)
grepl(longest.Dup4.end.sequence, right.Dup4)
# There is uncertainty about the exact location of the breakpoint. This is because there is a 10bp sequence that is 
# identical at the start point and the end point in the reference genome. Therefore, the reads align to it no matter
# which side of the breakpoint the aligner decides to put it. The breakpoint could be anywhere in that sequence. 
# Surely this can't be a coincidence. 
#
# GTGCACTCGCGGGAACTCGGACCTT   TCCATCGGGA   AGAACTCCTCCATCGTCGCAATCGTTG
#                         ^   ^        ^   ^
#         Position 28599070   ^        ^   Position 28595173
#  Position 28599071 / 28595163        Position 28595172 / 28599080

# Dup5
assumed.Dup5.start.position <- 28593122
assumed.Dup5.end.position <- 28598971
all.Dup5.start.sequences <- character()
for (sample.name in sub('-', '_', Dup5.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup5.start.position, 'Clipped_sequence']
	all.Dup5.start.sequences <- c(all.Dup5.start.sequences, these.start.sequences)
}
longest.Dup5.start.sequence <- all.Dup5.start.sequences[which.max(nchar(all.Dup5.start.sequences))]
# Do the same at the end of the sequence
all.Dup5.end.sequences <- character()
for (sample.name in sub('-', '_', Dup5.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup5.end.position, 'Clipped_sequence']
	all.Dup5.end.sequences <- c(all.Dup5.end.sequences, these.end.sequences)
}
longest.Dup5.end.sequence <- all.Dup5.end.sequences[which.max(nchar(all.Dup5.end.sequences))]
# There are no consistent alingment points for the reads at the start or end point. 
# Dup5 also includes a deletion from 28593930 to 28596533. This is supported by FM reads and BP reads. The BP reads 
# have a AGG that doesn't align at either end. 

# Dup5 deletion
assumed.Dup5del.start.position <- 28593930
assumed.Dup5del.end.position <- 28596533
all.Dup5del.start.sequences <- character()
for (sample.name in sub('-', '_', Dup5.samples)){
	these.start.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup5del.start.position, 'Clipped_sequence']
	all.Dup5del.start.sequences <- c(all.Dup5del.start.sequences, these.start.sequences)
}
longest.Dup5del.start.sequence <- all.Dup5del.start.sequences[which.max(nchar(all.Dup5del.start.sequences))]
# Do the same at the end of the sequence
all.Dup5del.end.sequences <- character()
for (sample.name in sub('-', '_', Dup5.samples)){
	these.end.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup5del.end.position, 'Clipped_sequence']
	all.Dup5del.end.sequences <- c(all.Dup5del.end.sequences, these.end.sequences)
}
longest.Dup5del.end.sequence <- all.Dup5del.end.sequences[which.max(nchar(all.Dup5del.end.sequences))]
# The deletion breakpoint looks like this
# CGTCTTCAAGGGTCAGCA         AGG         GTTGCACCACATCCGACGGA
#     start of del ^    inserted seq     ^ end of del
# position  28593929                     position 28596534

# Dup6: 
# There is only one Dup6 sample, but let's inspect it anyway
assumed.Dup6.start.position <- 28596241
assumed.Dup6.end.position <- 28602177
all.Dup6.start.sequences <- character()
for (sample.name in sub('-', '_', Dup6.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup6.start.position, 'Clipped_sequence']
	all.Dup6.start.sequences <- c(all.Dup6.start.sequences, these.start.sequences)
}
longest.Dup6.start.sequence <- all.Dup6.start.sequences[which.max(nchar(all.Dup6.start.sequences))]
# Do the same at the end of the sequence
all.Dup6.end.sequences <- character()
for (sample.name in sub('-', '_', Dup6.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup6.end.position, 'Clipped_sequence']
	all.Dup6.end.sequences <- c(all.Dup6.end.sequences, these.end.sequences)
}
longest.Dup6.end.sequence <- all.Dup6.end.sequences[which.max(nchar(all.Dup6.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point.
left.Dup6 <- substr(genome[4], assumed.Dup6.end.position-100, assumed.Dup6.end.position-1)
right.Dup6 <- substr(genome[4], assumed.Dup6.start.position+1, assumed.Dup6.start.position+100)
# The sequences are offset by one, this is because there is a C either side of the breakpoint, so we can't know
# for sure where the break happened, and the C aligns in both directions (ie: whichever side of the C is after the
# break, it still aligns because it's still a C
# We check that the clipped sequences are indeed as we expect. 
grepl(longest.Dup6.start.sequence, left.Dup6)
grepl(longest.Dup6.end.sequence, right.Dup6)
# Those two tests come out as false, because of a few SNPs, but I'm pretty sure we have the right sequences. 
# The breakpoint looks like this (the C is present at both ends and therefore could be either side of the breakpoint)
# GGCGGGTACTGTACAACA C GAAGCAATGCTGGCGATG
#   end of the dup ^   ^ start of the dup
#                  ^   ^ 
#  position 28602175   position 28596243

# Dup7:
assumed.Dup7.start.position <- 28597504
assumed.Dup7.end.position <- 28604250
all.Dup7.start.sequences <- character()
for (sample.name in sub('-', '_', Dup7.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup7.start.position, 'Clipped_sequence']
	all.Dup7.start.sequences <- c(all.Dup7.start.sequences, these.start.sequences)
}
longest.Dup7.start.sequence <- all.Dup7.start.sequences[which.max(nchar(all.Dup7.start.sequences))]
all.Dup7.end.sequences <- character()
for (sample.name in sub('-', '_', Dup7.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup7.end.position, 'Clipped_sequence']
	all.Dup7.end.sequences <- c(all.Dup7.end.sequences, these.end.sequences)
}
longest.Dup7.end.sequence <- all.Dup7.end.sequences[which.max(nchar(all.Dup7.end.sequences))]
# The longest.Dup7.end.sequence aligns to the same area of the unknown chromosome as the mates of the crosschrom reads. 
# For the start reads, the clipped bases before 28597504 can align to various different places in the genome (eg: 
# 2L:4665775, 2L:9360731, 2R:46566398). These region are all identical or nearly, but they only start being identical
# at the precise place where the clipped reads start aligning. 

#Dup8:
assumed.Dup8.start.position <- 28594797
assumed.Dup8.end.position <- 28602349
all.Dup8.start.sequences <- character()
for (sample.name in sub('-', '_', Dup8.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup8.start.position, 'Clipped_sequence']
	all.Dup8.start.sequences <- c(all.Dup8.start.sequences, these.start.sequences)
}
longest.Dup8.start.sequence <- all.Dup8.start.sequences[which.max(nchar(all.Dup8.start.sequences))]
# Do the same at the end of the sequence
all.Dup8.end.sequences <- character()
for (sample.name in sub('-', '_', Dup8.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup8.end.position, 'Clipped_sequence']
	all.Dup8.end.sequences <- c(all.Dup8.end.sequences, these.end.sequences)
}
longest.Dup8.end.sequence <- all.Dup8.end.sequences[which.max(nchar(all.Dup8.end.sequences))]
# If you look at the soft-clipped bases for both the start and end reads (did this in sample AC0159), they are actually 
# the SAME. This is probably why you get the crosschrom pairs matching to the same region on both sides of the 
# duplication. We find the same thing with dup5. 

# Dup9:
assumed.Dup9.start.position <- 28591140
assumed.Dup9.end.position <- 28601188
all.Dup9.start.sequences <- character()
for (sample.name in sub('-', '_', Dup9.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup9.start.position, 'Clipped_sequence']
	all.Dup9.start.sequences <- c(all.Dup9.start.sequences, these.start.sequences)
}
longest.Dup9.start.sequence <- all.Dup9.start.sequences[which.max(nchar(all.Dup9.start.sequences))]
# Do the same at the end of the sequence
all.Dup9.end.sequences <- character()
for (sample.name in sub('-', '_', Dup9.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup9.end.position, 'Clipped_sequence']
	all.Dup9.end.sequences <- c(all.Dup9.end.sequences, these.end.sequences)
}
longest.Dup9.end.sequence <- all.Dup9.end.sequences[which.max(nchar(all.Dup9.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point.
left.Dup9 <- substr(genome[4], assumed.Dup9.end.position-100, assumed.Dup9.end.position-1)
right.Dup9 <- substr(genome[4], assumed.Dup9.start.position+1, assumed.Dup9.start.position+100)
# We check that the clipped sequences are indeed as we expect. 
grepl(longest.Dup9.start.sequence, left.Dup9)
grepl(longest.Dup9.end.sequence, right.Dup9)
# Those two tests come out as false, but only because of there is a region of c.10 bases after / before the breakpoint
# where the reads don't align to any part of the expected genome sequence. Effectively, the breakpoint looks like this
# GCTGGACGAGTCCAAGTT  GATGAAGAAGAGAG  ATGAAATGTGTGCGTCATG
#   end of the dup ^    unkown seq    ^ start of the dup
#                  ^                  ^
#      position 28601187        position 28591144

# Dup10:
# Only two individuals appear to have Dup10, so it's hard to come to any firm conclusions. The first individual ('AN0099_C')
# has two FA read pairs supporting it, and two end breakpoint reads (pos 28603786). The second individual ('AN0047_C') has 
# three start breakpoint reads supporting it (pos 28593642). The clipped bases of these reads struggle to align anywhere 
# because a string of 27 bases (AAGACGCAGCGAATGTACTTTTTCGCT) sit between the end point (pos 28603786) and the subsequent 
# start point. AN0047_C has only start breakpoint reads, but AN0099 has end ones and they have the same inserted sequence, 
# confirming what is going on (and confirming that they are the same duplication). 
assumed.Dup10.start.position <- 28593642
assumed.Dup10.end.position <- 28603786
all.Dup10.start.sequences <- character()
for (sample.name in sub('-', '_', Dup10.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup10.start.position, 'Clipped_sequence']
	all.Dup10.start.sequences <- c(all.Dup10.start.sequences, these.start.sequences)
}
longest.Dup10.start.sequence <- all.Dup10.start.sequences[which.max(nchar(all.Dup10.start.sequences))]
# Do the same at the end of the sequence
all.Dup10.end.sequences <- character()
for (sample.name in sub('-', '_', Dup10.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup10.end.position, 'Clipped_sequence']
	all.Dup10.end.sequences <- c(all.Dup10.end.sequences, these.end.sequences)
}
longest.Dup10.end.sequence <- all.Dup10.end.sequences[which.max(nchar(all.Dup10.end.sequences))]
# Here is what the breakpoint looks like
# ACGCTTTTAATTTTATTC AAGACGCAGCGAATGTACTTTTTCGCT AAGAGGACGCAGCGGGTC\n',
#   end of the dup ^         inserted seq        ^ start of the dup\n',
#                  ^                             ^
#  position 28603785                             position 28593643\n',

# Dup11: 
assumed.Dup11.start.position <- 28581256
assumed.Dup11.end.position <- 28604994
all.Dup11.start.sequences <- character()
for (sample.name in sub('-', '_', Dup11.samples)){
	these.start.sequences <- clipping.end.point.gst.list[[sample.name]][clipping.end.point.gst.list[[sample.name]]$Position == assumed.Dup11.start.position, 'Clipped_sequence']
	all.Dup11.start.sequences <- c(all.Dup11.start.sequences, these.start.sequences)
}
longest.Dup11.start.sequence <- all.Dup11.start.sequences[which.max(nchar(all.Dup11.start.sequences))]
# Do the same at the end of the sequence
all.Dup11.end.sequences <- character()
for (sample.name in sub('-', '_', Dup11.samples)){
	these.end.sequences <- clipping.start.point.gst.list[[sample.name]][clipping.start.point.gst.list[[sample.name]]$Position == assumed.Dup11.end.position, 'Clipped_sequence']
	all.Dup11.end.sequences <- c(all.Dup11.end.sequences, these.end.sequences)
}
longest.Dup11.end.sequence <- all.Dup11.end.sequences[which.max(nchar(all.Dup11.end.sequences))]
# The segment upstream of the breakpoint (ie: the clipped sequences on the left of the reads that align at the start
# of the duplicated segment) is present in more than one place on the unknown chromosome (before 29210988, 32073377) 
# and also on other chromosomes (before X:24393108, before 2L:4376816). These are just three places I found, there 
# may be more. Since in most samples with this duplication there are crosschrom pairs on both sides of the duplication 
# that align to the same point on UNKN (~29210980). 

##################
# Genotype calls #
##################

# Create a function that will take a logical vector and output the indices of the ranges of consecutive True (T)
# values, removing a set number of indices either side as a buffer. 
T.windows <- function(x, window.buffer = 0){
	T.so.far <- 0
	output <- matrix(0,0,3)
	colnames(output) <- c('start', 'end', 'length')
	start.pos <- 0
	for (i in 1:length(x)){
		if (x[i]){
			if (start.pos == 0){
				if (((T.so.far == window.buffer) | (i ==1))) # If this is the first position of the vector, we don't need the buffer at the beginning
					start.pos <- i
				else
					T.so.far <- T.so.far + 1
			}
		}
		else {
			if (start.pos){
				this.output <- c(start.pos, i - 1 - window.buffer)
				length.of.T <- this.output[2] - this.output[1] + 1
				if (length.of.T > 0){
					output <- rbind(output, c(this.output, length.of.T))
				}
				start.pos <- 0
			}
			T.so.far <- 0
		}
	}
	# For the last window, we don't need the buffer on the right hand side
	if (start.pos){
		this.output <- c(start.pos, i)
		length.of.T <- this.output[2] - this.output[1] + 1
		if (length.of.T > 0)
			output <- rbind(output, c(this.output, length.of.T))
	}
	return(output)
}


# Write a function that makes a genotype call for a given duplication. "window.buffer" is the number of windows 
# we trim off the edge of a range of windows to calculate mean coverage. "minimum.range" is the minimum number 
# of windows, after trimming, that can be used to estimate coverage
make.genotype.calls <- function(sample.name, method = 'new', window.buffer = 0, minimum.windows = 3, verbose = F){
	if (! method %in% c('orig', 'new'))
		stop('"method" should be one of "orig" or "new".')
	if (method == 'new')
		counts.table <- counts.list[[sample.name]]
	else
		counts.table <- counts.list_orig[[sample.name]]
	sample.name <- sub('-', '_', sample.name)
	if (verbose) cat('Sample: ', sample.name, '\n')
	# This is the vector that we are going to return after modification
	return.vector <- rep(2,ncol(read.based.gst.duplications))
	names(return.vector) <- colnames(read.based.gst.duplications)
	# Manually determine genomic regions for each of the duplication types:
	Dup1.region <- Dup1.bp
	Dup2.region <- Dup2.bp
	# The region of Dup3 doesn't stretch as far as the breakpoints suggest. 
	Dup3.region <- c(Dup3.bp[1], 28592100)
	Dup4.region <- Dup4.bp
	# There is a deletion in Dup5, so restrict the region to where the deletion isn't
	Dup5.region <- c(28596300, Dup5.bp[2])
	Dup6.region <- Dup6.bp
	Dup7.region <- Dup7.bp
	Dup8.region <- Dup8.bp
	Dup9.region <- Dup9.bp
	Dup10.region <- Dup10.bp
	Dup11.region <- Dup11.bp
	# Create regions as logical vectors for each of those
	pos.list <- list()
	pos.list[['Dup1']] <- (counts.table$Position >= Dup1.region[1]) & (counts.table$Position <= Dup1.region[2] - coverage.window)
	pos.list[['Dup2']] <- (counts.table$Position >= Dup2.region[1]) & (counts.table$Position <= Dup2.region[2] - coverage.window)
	pos.list[['Dup3']] <- (counts.table$Position >= Dup3.region[1]) & (counts.table$Position <= Dup3.region[2] - coverage.window)
	pos.list[['Dup4']] <- (counts.table$Position >= Dup4.region[1]) & (counts.table$Position <= Dup4.region[2] - coverage.window)
	pos.list[['Dup5']] <- (counts.table$Position >= Dup5.region[1]) & (counts.table$Position <= Dup5.region[2] - coverage.window)
	pos.list[['Dup6']] <- (counts.table$Position >= Dup6.region[1]) & (counts.table$Position <= Dup6.region[2] - coverage.window)
	pos.list[['Dup7']] <- (counts.table$Position >= Dup7.region[1]) & (counts.table$Position <= Dup7.region[2] - coverage.window)
	pos.list[['Dup8']] <- (counts.table$Position >= Dup8.region[1]) & (counts.table$Position <= Dup8.region[2] - coverage.window)
	pos.list[['Dup9']] <- (counts.table$Position >= Dup9.region[1]) & (counts.table$Position <= Dup9.region[2] - coverage.window)
	pos.list[['Dup10']] <- (counts.table$Position >= Dup10.region[1]) & (counts.table$Position <= Dup10.region[2] - coverage.window)
	pos.list[['Dup11']] <- (counts.table$Position >= Dup11.region[1]) & (counts.table$Position <= Dup11.region[2] - coverage.window)
	#
	# We write a function to identify the positions unique to each of the duplication and calculate coverage 
	# at those positions
	find.cnv.state.from.unique.positions <- function(rd, ad, rv, v){
		unique.pos <- list()
		has.unique.windows <- character()
		output.rv <- rv
		# Find dups in this sample that have unique windows
		for	(this.dup in rd){
			this.dup.pos <- pos.list[[this.dup]]
			other.dups <- rd[rd != this.dup]
			for (o in other.dups)
				this.dup.pos <- this.dup.pos & (!pos.list[[o]])
			# Get the ranges of consective True values in the vector 
			unique.pos[[this.dup]] <- T.windows(this.dup.pos, window.buffer)
			if (sum(unique.pos[[this.dup]][,'length']) > minimum.windows)
				has.unique.windows <- c(has.unique.windows, this.dup)
		}
		# If none of the dups have unique windows, we can't make the call
		if (length(has.unique.windows) == 0){
			output.rv[rd] <- NA
			if (v) cat('\t\tFailed to make calls on sample ', sample.name, ' for ', paste(rd), ' because no duplications had unique windows.\n', sep = '')
			return(list(result = F, rv = output.rv))
		}
		# Otherwise, make the calls for the dups with unique windows 
		for (this.dup in has.unique.windows){
			all.unique.pos.vector <- numeric()
			u <- unique.pos[[this.dup]]
			for (i in 1:nrow(u))
				all.unique.pos.vector <- c(all.unique.pos.vector, u[i,1]:u[i,2])
			# Now for each position of this dup, subtract the coverage due to the coverage at other dups
			this.dup.temp.coverage <- counts.table[,'CNV']
			call.failed <- F
			for (this.other.dup in ad){
				# Check whether these two dups overlap
				if (any((pos.list[[this.dup]] + pos.list[[this.other.dup]]) == 2)){
					# Check whether the coverage for the overlapping dup is NA
					if (is.na(return.vector[this.other.dup])){
						output.rv[this.dup] <- NA
						if (v) cat('\t\tFailed to make call on sample ', sample.name, ' for ', this.dup, ' because it overlaps with ', this.other.dup, ', which has NA coverage.\n', sep = '')
						call.failed <- T
						break
					}
					# Otherwise, subtract the coverage of the overlapping dup
					this.dup.temp.coverage <- this.dup.temp.coverage - pos.list[[this.other.dup]]*(output.rv[this.other.dup] - 2)
				}
			}
			if (call.failed)
				next
			# Get the most frequent CNV call for the positions of interest. 
			cov.table <- sort(table(this.dup.temp.coverage[all.unique.pos.vector]), decreasing = T)
			# Make the call if that most frequent coverage has at least 70% frequency
			if ((cov.table[1]/sum(cov.table)) >= 0.7)
				median.cov <- as.numeric(names(cov.table)[1])
			else{
				if (v) cat('\t\tFailed to make call on sample ', sample.name, ' for ', this.dup, ' because coverage was too variable.\n', sep = '')
				median.cov <- NA
			}
			output.rv[this.dup] <- median.cov
		}
		return(list(result = T, rv = output.rv, has.unique.windows = has.unique.windows))
	}
	#
	# OK, now we check which duplications the sample carries
	sample.dups <- colnames(read.based.gst.duplications)[read.based.gst.duplications[sample.name, ] >= 1]
	sample.dups <- sample.dups[sample.dups != 'Dup0']
	remaining.dups <- sample.dups
	analysed.dups <- character()
	# If we have a Dup0, then there is something going on that we don't have a handle on so we bail
	Dup0 <- read.based.gst.duplications[sample.name, 'Dup0']
	if (Dup0 > 0){
		if (verbose){
			if (length(sample.dups))
				cat('\tFailed to make calls on sample ', sample.name, ' because there was an unidentified duplication.\n', sep = '')
			else
				cat('\t', sample.name, ' had an unidentified duplication, but no others.\n', sep = '')
		}
		return.vector[sample.dups] <- NA
		return(return.vector)
	}
	if (verbose & (length(sample.dups) == 0)) cat('\tNo duplications detected in this sample\n')
	# Keep running the main function until there are no dups left with unique windows
	while(length(remaining.dups)){
		if (verbose) cat('\tThe following dups remain to be processed: ', paste(remaining.dups, collapse = ';'), '\n', sep = '')
		unique.pos.cnv <- find.cnv.state.from.unique.positions(remaining.dups, analysed.dups, return.vector, verbose)
		return.vector <- unique.pos.cnv[['rv']]
		if (unique.pos.cnv[['result']]){
			analysed.dups <- c(analysed.dups, unique.pos.cnv[['has.unique.windows']])
			remaining.dups <- setdiff(sample.dups, analysed.dups)
		}
		else 
			return(return.vector)
	}
	return(return.vector)
}

# Get the genotypes for all the samples
genotype.call.1 <- make.genotype.calls(sample.names[1], verbose = T)
genotype.calls <- matrix(0, nrow = length(sample.names), ncol = length(genotype.call.1), dimnames = list(sample.names, names(genotype.call.1)))
genotype.calls[sample.names[1], ] <- genotype.call.1
for (s in sample.names[2:length(sample.names)])
	genotype.calls[s, ] <- make.genotype.calls(s, verbose = T)

# Let's have a look at the distribution of Dup1:
Dup1.pop.table <- table(genotype.calls[,'Dup1'], meta$population)
Dup1.BFM.reg.table <- table(genotype.calls[meta$population == 'BFcol','Dup1'], as.factor(as.character(meta$region[meta$population == 'BFcol'])))
HWExact(as.numeric(Dup1.pop.table[,'BFcol'])) 

# Dup2:
Dup2.pop.table <- table(genotype.calls[,'Dup2'], meta$population)
Dup2.UG.reg.table <- table(genotype.calls[meta$population == 'UGgam','Dup2'], as.factor(as.character(meta$region[meta$population == 'UGgam'])))

# Dup3:
Dup3.pop.table <- table(genotype.calls[,'Dup3'], meta$population)
HWExact(as.numeric(Dup3.pop.table[,'UGgam'])) 

# Dup4:
Dup4.pop.table <- table(genotype.calls[,'Dup4'], meta$population)

# Dup5:
Dup5.pop.table <- table(genotype.calls[,'Dup5'], meta$population)

# Dup6:
Dup6.pop.table <- table(genotype.calls[,'Dup6'], meta$population)

# Dup7:
Dup7.pop.table <- table(genotype.calls[,'Dup7'], meta$population)

# Dup8:
Dup8.pop.table <- table(genotype.calls[,'Dup8'], meta$population)

# Dup9:
Dup9.pop.table <- table(genotype.calls[,'Dup9'], meta$population)
HWExact(c(as.numeric(Dup9.pop.table[,'AOcol']), 0)) 

# Dup10:
Dup10.pop.table <- table(genotype.calls[,'Dup10'], meta$population)

# Dup11:
Dup11.pop.table <- table(genotype.calls[,'Dup11'], meta$population)
HWExact(as.numeric(Dup11.pop.table[,'KE'])) 

write.table(genotype.calls, 'genotype_calls.csv', sep = '\t', col.names = NA)

save.image('Gste2_analysis_shrunk_data.Rdata')



