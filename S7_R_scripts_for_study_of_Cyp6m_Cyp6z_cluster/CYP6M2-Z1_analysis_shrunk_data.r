library(stringr)
library(Biostrings)
library(HardyWeinberg)

load('CNV_analysis_phase2_for_CYP6M2-Z1.Rdata')

# Get the coordinates of all the genes in the cluster
cyp6m2.coords <- c(6928858, 6928858+1628)
cyp6m3.coords <- c(6931811, 6931811+1888)
cyp6m4.coords <- c(6933589, 6933589+1862)
cyp6z1.coords <- c(6976501, 6976501+1641)
cyp6z2.coords <- c(6973779, 6973779+1784)
cyp6z3.coords <- c(6971669, 6971669+1621)

cyp6m2.z1.region.coord <- c(6900000, 7000000) 
cyp6m2.z1.region.plotting.indices <- c(min(which(counts.list[[1]]$Position >= cyp6m2.z1.region.coord[1])), max(which(counts.list[[1]]$Position <= cyp6m2.z1.region.coord[2])))

# Identify the samples with a detected duplication in cyp6m2 and cyp6z1
encompass.duplication.detected.in.cyp6m2 <- encompass.output.by.gene[['AGAP008212']]
encompass.duplication.detected.in.cyp6z1 <- encompass.output.by.gene[['AGAP008219']]
encompass.duplication.detected.in.either <- unique(c(encompass.duplication.detected.in.cyp6m2, 
                                                         encompass.duplication.detected.in.cyp6z1))

# Load table containing geographical id
meta <- read.table('/home/eric/Liverpool/phase2.AR1/samples/samples.meta.txt', header = T, row.names = 1, sep = '\t', quote = "", comment.char = "")
rownames(meta) <- sub('-', '_', rownames(meta))
# Create booleans vector that tell us whether each individual in the samples table has a duplication in cyp6aa1
# or cyp6p3
hasencomp.cyp6m2 <- rownames(meta) %in% encompass.duplication.detected.in.cyp6m2
hasencomp.cyp6z1 <- rownames(meta) %in% encompass.duplication.detected.in.cyp6z1
# add these data to the meta table
meta$hasencomp.cyp6m2 = hasencomp.cyp6m2
meta$hasencomp.cyp6z1 = hasencomp.cyp6z1

# Now, within each population, we want to know the proportion of individuals carrying each of those duplications
propencomp.cyp6m2 <- tapply(meta$hasencomp.cyp6m2, meta$population, mean)
propencomp.cyp6z1 <- tapply(meta$hasencomp.cyp6z1, meta$population, mean)
propdup.for.plot.cyp6m2 <- rbind(propencomp.cyp6m2, 1 - propencomp.cyp6m2)
propdup.for.plot.cyp6z1 <- rbind(propencomp.cyp6z1, 1 - propencomp.cyp6z1)
pop.sizes <- tapply(meta$population, meta$population, length)

x11(width = 15)
par(mfrow = c(2,1), mar = c(2,4,4,8), oma = c(3,0,0,0), xpd = F)
barplot(propdup.for.plot.cyp6m2, legend.text = c('Duplication', 'No duplication'), ylab = 'Proportion', names.arg = rep('', ncol(propdup.for.plot.cyp6m2)), main = 'CYP6M2', args.legend = list(x=22.2, y=0.9))
barplot(propdup.for.plot.cyp6z1, ylab = 'Proportion', names.arg = rep('', ncol(propdup.for.plot.cyp6z1)),  main = 'CYP6Z1')
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
SSFA.folder <- paste('/home/eric/Liverpool/CNV_/SSFA/v3_', chrom, '/phase2/CYP6M2-Z1_region', sep = '')
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

# Get the faceaway reads in the region
FA.cyp.list <- list()
for (i in 1:length(FA.list)){
	this.list <- FA.list[[i]]
	this.name <- names(FA.list)[i]
	# Use the following line if you want reads that start OR end in the region
	which.in.cyp.FA <- ((this.list$Position >= cyp6m2.z1.region.coord[1]) & (this.list$Position <= cyp6m2.z1.region.coord[2])) | ((this.list$Mate.position >= cyp6m2.z1.region.coord[1]) & (this.list$Mate.position <= cyp6m2.z1.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the region
#	which.in.cyp.FA <- ((this.list$Position >= cyp6m2.z1.region.coord[1]) & (this.list$Mate.position <= cyp6m2.z1.region.coord[2]))
	FA.cyp.list[[this.name]] <- this.list[which.in.cyp.FA, ]
}
SS.cyp.list <- list()
for (i in 1:length(SS.list)){
	this.list <- SS.list[[i]]
	this.name <- names(SS.list)[i]
	# Use the following line if you want reads that start OR end in the region
	which.in.cyp.SS <- ((this.list$Position >= cyp6m2.z1.region.coord[1]) & (this.list$Position <= cyp6m2.z1.region.coord[2])) | ((this.list$Mate.position >= cyp6m2.z1.region.coord[1]) & (this.list$Mate.position <= cyp6m2.z1.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the region
#	which.in.cyp.SS <- ((this.list$Position >= cyp6m2.z1.region.coord[1]) & (this.list$Mate.position <= cyp6m2.z1.region.coord[2]))
	SS.cyp.list[[this.name]] <- this.list[which.in.cyp.SS, ]
}
FM.cyp.list <- list()
for (i in 1:length(FM.list)){
	this.list <- FM.list[[i]]
	this.name <- names(FM.list)[i]
	# Use the following line if you want reads that start OR end in the region
	which.in.cyp.FM <- ((this.list$Position >= cyp6m2.z1.region.coord[1]) & (this.list$Position <= cyp6m2.z1.region.coord[2])) | ((this.list$Mate.position >= cyp6m2.z1.region.coord[1]) & (this.list$Mate.position <= cyp6m2.z1.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the region
#	which.in.cyp.FM <- ((this.list$Position >= cyp6m2.z1.region.coord[1]) & (this.list$Mate.position <= cyp6m2.z1.region.coord[2]))
	FM.cyp.list[[this.name]] <- this.list[which.in.cyp.FM, ]
}
crosschrom.cyp.list <- list()
for (i in 1:length(crosschrom.list)){
	this.list <- crosschrom.list[[i]]
	this.name <- names(crosschrom.list)[i]
	which.in.cyp.crosschrom <- (this.list$Position >= cyp6m2.z1.region.coord[1]) & (this.list$Position <= cyp6m2.z1.region.coord[2])
	crosschrom.cyp.list[[this.name]] <- this.list[which.in.cyp.crosschrom, ]
}

# Load up the results of the breakpoint detection 
cat('Loading breakpoints data\n')
bp.folder <- paste('/home/eric/Liverpool/CNV_/breakpoints/' , chrom, '/phase2/CYP6M2-Z1_region', sep = '')
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
clipping.start.point.cyp.list <- list()
for (i in 1:length(clipping.start.point.list)){
	this.list <- clipping.start.point.list[[i]]
	this.name <- names(clipping.start.point.list)[i]
	which.in.cyp.clipping.start.point <- (this.list$Position >= cyp6m2.z1.region.coord[1]) & (this.list$Position <= cyp6m2.z1.region.coord[2])
	this.initial.list <- this.list[which.in.cyp.clipping.start.point, , drop = F]
	# The following for lines allow breakpoints to be kept only if their frequency exceeds a certain threshold. They 
	# are effectively pointless unless you set the threshold below above 0.
	reads.per.start.point <- tapply(this.initial.list$Position, as.factor(this.initial.list$Position), length)
	positions.to.keep <- this.initial.list$Position %in% names(reads.per.start.point)[reads.per.start.point > 0]
	clipping.start.point.cyp.list[[this.name]] <- this.initial.list[positions.to.keep, ]
}
clipping.end.point.cyp.list <- list()
for (i in 1:length(clipping.end.point.list)){
	this.list <- clipping.end.point.list[[i]]
	this.name <- names(clipping.end.point.list)[i]
	which.in.cyp.clipping.end.point <- (this.list$Position >= cyp6m2.z1.region.coord[1]) & (this.list$Position <= cyp6m2.z1.region.coord[2])
	this.initial.list <- this.list[which.in.cyp.clipping.end.point, , drop = F]
	reads.per.end.point <- tapply(this.initial.list$Position, as.factor(this.initial.list$Position), length)
	positions.to.keep <- this.initial.list$Position %in% names(reads.per.end.point)[reads.per.end.point > 0]
	clipping.end.point.cyp.list[[this.name]] <- this.initial.list[positions.to.keep, ]
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
mapq0.proportions <- subset(mapq0.proportions, (Position >= cyp6m2.z1.region.coord[1]) & (Position <= cyp6m2.z1.region.coord[2]))
total.coverage <- mapq0.proportions$'Count mapq > 0' + mapq0.proportions$'Count mapq = 0'
mapq0.proportions$mapq0.prop <- mapq0.proportions$'Count mapq = 0' / total.coverage
# In some cases, there is no coverage, so we get infinite values. Let's deal with those
mapq0.proportions$mapq0.prop[total.coverage == 0] <- 0

# Create a list of colours that will be used for plotting the duplication types
duplication.colours <- c(rgb(0.6509804, 0.8078431, 0.8901961), rgb(0.1215686, 0.4705882, 0.7058824))#, rgb(0.6980392, 0.8745098, 0.5411765), rgb(0.2, 0.627451, 0.172549), rgb(0.9843137, 0.6039216, 0.6), rgb(0.8901961, 0.1019608, 0.1098039), rgb(0.9921569, 0.7490196, 0.4352941), rgb(1, 0.4980392, 0), rgb(0.7921569, 0.6980392, 0.8392157), rgb(0.4156863, 0.2392157, 0.6039216), rgb(1,0,1,0.7), rgb(0.6941176, 0.3490196, 0.1568627), rgb(0.5,0.5,0,0.8), rgb(0,0.5,0.5,0.8), rgb(0.5, 0, 0.5, 0.8))

plot.cyp <- function(this.sample, list.of.dups = NULL, FA = T, SS=T, FM=T, XC=T, BP=T, start.index = cyp6m2.z1.region.plotting.indices[1], end.index = cyp6m2.z1.region.plotting.indices[2]){
	this.sample.name <- regmatches(this.sample, regexpr('A.\\d\\d\\d\\d_C', this.sample))
	plot(counts.list[[this.sample]]$Position[start.index : end.index], counts.list[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		if (list.of.dups['Dupz1'] >= 1){
			rect(6968962, -0.6, 6979681, 0, col = rgb(1, 0, 0, 0.8), border = rgb(1, 0, 0, 0.8))
			text(6979681, -0.4, 'Dupz1', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dupm2'] >= 1){
			rect(6927942, -0.6, 6930500, 0, col = rgb(0, 1, 0, 0.8), border = rgb(0, 1, 0, 0.8))
			text(6930500, -0.4, 'Dupm2', cex = 1.2, adj = c(1,0))
		}
	}
	if (FA){
		if (nrow(FA.cyp.list[[this.sample.name]]) > 0)
			add.faceaways(FA.cyp.list[[this.sample.name]], this.col = 'blue', yrange = c(0,3))
	}
	if (SS){
		if (nrow(SS.cyp.list[[this.sample.name]]) > 0)
			add.faceaways(SS.cyp.list[[this.sample.name]], this.col = 'cyan', yrange = c(3,6))
	}
	if (FM){
		if (nrow(FM.cyp.list[[this.sample.name]]) > 0)
			add.faceaways(FM.cyp.list[[this.sample.name]], this.col = 'green', yrange = c(6,9))
	}
	if (XC){
		if (nrow(crosschrom.cyp.list[[this.sample.name]])) 
			add.breakpoints(crosschrom.cyp.list[[this.sample.name]]$Position, this.col = 'red', yrange = c(9,12))
	}
	if (BP){
		if (nrow(clipping.end.point.cyp.list[[this.sample.name]]))
			add.breakpoints(clipping.end.point.cyp.list[[this.sample.name]][, 'Position'], this.col = 'brown', yrange = c(0,12))
		if (nrow(clipping.start.point.cyp.list[[this.sample.name]]))
			add.breakpoints(clipping.start.point.cyp.list[[this.sample.name]][, 'Position'], this.col = 'pink', yrange = c(0,12))
	}
	these.cnv.states <- counts.list[[this.sample]]$CNV[start.index : end.index]
	lines(counts.list[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2)
	abline(v = cyp6m2.coords[1])
	abline(v = cyp6m2.coords[2])
	text(mean(cyp6m2.coords), 9, 'CYP6M2', srt=90, adj = 0)
	# Let's show other genes as well
	abline(v = cyp6m3.coords[1], col = 'brown')
	abline(v = cyp6m3.coords[2], col = 'brown')
	text(mean(cyp6m3.coords), 9, 'CYP6M3', srt=90, adj = 0, col = 'brown')
	abline(v = cyp6m4.coords[1], col = 'violet')
	abline(v = cyp6m4.coords[2], col = 'violet')
	text(mean(cyp6m4.coords), 9, 'CYP6M4', srt=90, adj = 0, col = 'violet')
	abline(v = cyp6z1.coords[1], col = 'green')
	abline(v = cyp6z1.coords[2], col = 'green')
	text(mean(cyp6z1.coords), 9, 'CYP6Z1', srt=90, adj = 0, col = 'green')
	abline(v = cyp6z2.coords[1], col = 'blue')
	abline(v = cyp6z2.coords[2], col = 'blue')
	text(mean(cyp6z2.coords), 9, 'CYP6Z2', srt=90, adj = 0, col = 'blue')
	abline(v = cyp6z3.coords[1], col = 'red')
	abline(v = cyp6z3.coords[2], col = 'red')
	text(mean(cyp6z3.coords), 9, 'CYP6Z3', srt=90, adj = 0, col = 'red')
}


plot.all.cyp <- function(list.of.all, list.of.dupl = encompass.duplication.detected.in.either, FA=T, SS=T, FM=T, XC=T, BP=T, start.index = cyp6m2.z1.region.plotting.indices[1], end.index = cyp6m2.z1.region.plotting.indices[2]){
	x11()
	list.of.all <- sub('-', '_', list.of.all)
	list.of.dupl <- sub('-', '_', list.of.dupl)
	for (this.sample in list.of.all){
		if (this.sample %in% list.of.dupl)
			par(bg = 'beige')
		else
			par(bg = 'white')
		plot.cyp(this.sample, NULL, FA, SS, FM, XC, BP, start.index, end.index)
		locator(1)
	}
}

plot.all.cyp.with.read.dups <- function(list.of.samples, matrix.of.read.dups = read.based.cyp.duplications, FA = T, SS = T, FM = T, XC=T, BP=T, start.index = cyp6m2.z1.region.plotting.indices[1], end.index = cyp6m2.z1.region.plotting.indices[2]){
	x11()
	x.midpoint <- mean(c(counts.list[[1]]$Position[start.index], counts.list[[1]]$Position[end.index]))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- sub('-', '_', list.of.samples[i])
		plot.cyp(this.sample, matrix.of.read.dups[this.sample,], FA, SS, FM, XC, BP, start.index, end.index)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


# Use the following definitions for the different duplications (numbers represent start range (1st row) and end range
# (2nd row) of the read pairs). These values were determined after visual inspection of the data. 
Dupz1.fa <- data.frame(c(6968950, 6979300), c(6969250, 6979600))
Dupz1.bp = c(6968962, 6979681)
Dupz1.bp.seq = c('ACGCT', 'AGGTT')
Dupm2.bp = c(6927942)
Dupm2.bp.seq = c('ATTAT')

# We will store the duplication type information in a matrix
read.based.cyp.duplications <- matrix(0, length(counts.list), 3, dimnames = list(names(counts.list), c('Dup0', 'Dupz1', 'Dupm2')))
# Now go through each individual and determine which duplications they have based on their diagnostic reads. 
for (this.sample in names(counts.list)){
	# get the FA and SS reads encompassed in the general region
	these.FA <- FA.cyp.list[[this.sample]]
	these.SS <- SS.cyp.list[[this.sample]]
	these.CSP <- clipping.start.point.cyp.list[[this.sample]]
	these.CEP <- clipping.end.point.cyp.list[[this.sample]]
	# Count the ones that start within the Dup1 start region and end within the Dup1 end region
	num.Dupz1.fa.reads <- sum((these.FA$Position > Dupz1.fa[1,1]) & these.FA$Position < Dupz1.fa[1,2]
	                       & (these.FA$Mate.position > Dupz1.fa[2,1]) & (these.FA$Mate.position < Dupz1.fa[2,2]))
	Dupz1.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dupz1.bp[1]]
	num.Dupz1.start.reads <- sum(substr(reverse(Dupz1.CEP.seq), 1, 5) == Dupz1.bp.seq[1])
	Dupz1.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dupz1.bp[2]]
	num.Dupz1.end.reads <- sum(substr(Dupz1.CSP.seq, 1, 5) == Dupz1.bp.seq[2])
	num.Dupz1.reads <- num.Dupz1.fa.reads + num.Dupz1.start.reads + num.Dupz1.end.reads
	#
	Dupm2.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dupm2.bp[1]]
	num.Dupm2.start.reads <- sum(substr(reverse(Dupm2.CEP.seq), 1, 5) == Dupm2.bp.seq[1])
	num.Dupm2.reads <- num.Dupm2.start.reads
	# For each of these, if there are at least 2 supporting reads, score it as a duplication. If there is only one
	# supporting read, score it tentatively
	read.based.cyp.duplications[this.sample, 2:3] <- as.logical(c(num.Dupz1.reads, num.Dupm2.reads))
	if (num.Dupz1.reads == 1)
		read.based.cyp.duplications[this.sample, 'Dupz1'] <- 0.5
	if (num.Dupm2.reads == 1)
		read.based.cyp.duplications[this.sample, 'Dupm2'] <- 0.5
	# If the sample has a duplication, but is none of the specific ones, we record it as Dup0
	if ((sum(read.based.cyp.duplications[this.sample, ] >= 1) == 0) & (this.sample %in% c(encompass.duplication.detected.in.either)))
		read.based.cyp.duplications[this.sample, 'Dup0'] <- 1
}

# Compared the duplications detected by coverage or reads
cov.based.duplications.names <- encompass.duplication.detected.in.either
read.based.duplications.names <- rownames(read.based.cyp.duplications)[as.logical(apply(read.based.cyp.duplications[,2:3] >= 1, 1, sum))]

cov.based.negatives <- setdiff(read.based.duplications.names, cov.based.duplications.names)
read.based.negatives <- setdiff(cov.based.duplications.names, read.based.duplications.names)
# We've now totally nailed the duplication calling with discordant reads. 

# Output a table of each duplication type and the country they are found in
Dup0.meta <- meta[rownames(read.based.cyp.duplications)[read.based.cyp.duplications[,'Dup0'] >= 1], 'population', drop = F]
Dupz1.meta <- meta[rownames(read.based.cyp.duplications)[read.based.cyp.duplications[,'Dupz1'] >= 1], 'population', drop = F]
Dupm2.meta <- meta[rownames(read.based.cyp.duplications)[read.based.cyp.duplications[,'Dupm2'] >= 1], 'population', drop = F]

# Let's output some lists of individuals carrying each duplication
write(rownames(Dup0.meta), file = 'CYP6_Dup0_samples.txt', ncolumns = 1)
write(rownames(Dupz1.meta), file = 'CYP6_Dupz1_samples.txt', ncolumns = 1)
write(rownames(Dupm2.meta), file = 'CYP6_Dupm2_samples.txt', ncolumns = 1)

Dup0.samples <- rownames(Dup0.meta)
Dupz1.samples <- rownames(Dupz1.meta)
Dupm2.samples <- rownames(Dupm2.meta)


######################################
# PRECISELY IDENTIFYING BREAK POINTS #
######################################

# Load the genome
genome <- readDNAStringSet('/home/eric/Liverpool/AR3/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa')

# Looking at the plot for individual genes, it wasn't immediately obvious that there were any breakpoint reads.
# I wonder if they would stand out more if I grouped all the samples together
clipping.end.point.allz1dups <- do.call(rbind, clipping.end.point.list[sub('-', '_', encompass.duplication.detected.in.cyp6z1)])
clipping.end.point.allz1dups <- clipping.end.point.allz1dups[order(clipping.end.point.allz1dups$Position), ]
clipping.start.point.allz1dups <- do.call(rbind, clipping.start.point.list[sub('-', '_', encompass.duplication.detected.in.cyp6z1)])
clipping.start.point.allz1dups <- clipping.start.point.allz1dups[order(clipping.start.point.allz1dups$Position), ]
# Plot these points on a graph
plot(c(counts.list[[1]]$Position[cyp6m2.z1.region.plotting.indices[1]], counts.list[[1]]$Position[cyp6m2.z1.region.plotting.indices[2]]), c(0, 0), ylim = c(0,6), xlim = c(6965000, 6985000), type = 'n')
add.breakpoints(clipping.end.point.allz1dups[, 'Position'], this.col = 'brown', yrange = c(0,6))
add.breakpoints(clipping.start.point.allz1dups[, 'Position'], this.col = 'pink', yrange = c(0,6))
# Add the rough points where the FA reads start and end
abline(v = Dupz1.fa[1,1])
abline(v = Dupz1.fa[2,2])
# That didn't help 

# Let's look more closely where the breakpoints are expected to be. 
clipping.end.point.allz1dups.narrow <- clipping.end.point.allz1dups[(clipping.end.point.allz1dups$Position > (Dupz1.fa[1,1] - 300)) & (clipping.end.point.allz1dups$Position < Dupz1.fa[1,2]), ]
plot(c(Dupz1.fa[1,1] - 300, Dupz1.fa[1,2]), c(0, 0), ylim = c(0,6), type = 'n')
add.breakpoints(clipping.end.point.allz1dups.narrow[, 'Position'], this.col = 'brown', yrange = c(0,6))
abline(v = Dupz1.fa[1,1])
abline(v = Dupz1.fa[1,2])
startpostable.z1 <- table(clipping.end.point.allz1dups.narrow$Position)
startpostable.z1 <- startpostable.z1[order(startpostable.z1, decreasing = T)]
# OK, we've got something for the start point! It's at position 6968962 
clipping.start.point.allz1dups.narrow <- clipping.start.point.allz1dups[(clipping.start.point.allz1dups$Position > Dupz1.fa[2,1]) & (clipping.start.point.allz1dups$Position < (Dupz1.fa[2,2] + 300)), ]
plot(c(Dupz1.fa[2,1], Dupz1.fa[2,2] + 300), c(0, 0), ylim = c(0,6), type = 'n')
add.breakpoints(clipping.start.point.allz1dups.narrow[, 'Position'], this.col = 'pink', yrange = c(0,6))
abline(v = Dupz1.fa[2,1])
abline(v = Dupz1.fa[2,2])
endpostable.z1 <- table(clipping.start.point.allz1dups.narrow$Position)
endpostable.z1 <- endpostable.z1[order(endpostable.z1, decreasing = T)]
# We have three candidates there: 6979535, 6979650 and 6979681
# Let's compare the frequency of these positions between individuals that do or don't have the dup (according
# to coverage).
# First the start point
z1 <- sub('-', '_', encompass.duplication.detected.in.cyp6z1)
nonz1 <- setdiff(names(counts.list), z1)
z1.scounts <- numeric()
nonz1.scounts <- numeric()
for (spos in names(startpostable.z1)){
	this.z1.count <- 0
	this.nonz1.count <- 0
	for (this.sample in z1)
		this.z1.count <- this.z1.count + any(clipping.end.point.list[[this.sample]]$Position == spos)
	for (this.sample in nonz1)
		this.nonz1.count <- this.nonz1.count + any(clipping.end.point.list[[this.sample]]$Position == spos)
	z1.scounts[as.character(spos)] <- this.z1.count
	nonz1.scounts[as.character(spos)] <- this.nonz1.count
}
all.scounts.z1 <- cbind(z1.scounts, nonz1.scounts)
all.scounts.z1 <- all.scounts.z1[order(all.scounts.z1[,1], decreasing = T),]
# Our best candidate (6968962) is also present in 6 apparent nonz1 samples. One of them (AB0178) is a false 
# negative. For the other 5, if you look at the clipped sequence, it is not the same as in our true positives, 
# so that's ok. 

# Now the end point
z1.ecounts <- numeric()
nonz1.ecounts <- numeric()
for (epos in names(endpostable.z1)){
	this.z1.count <- 0
	this.nonz1.count <- 0
	for (this.sample in z1)
		this.z1.count <- this.z1.count + any(clipping.start.point.list[[this.sample]]$Position == epos)
	for (this.sample in nonz1)
		this.nonz1.count <- this.nonz1.count + any(clipping.start.point.list[[this.sample]]$Position == epos)
	z1.ecounts[as.character(epos)] <- this.z1.count
	nonz1.ecounts[as.character(epos)] <- this.nonz1.count
}
all.ecounts.z1 <- cbind(z1.ecounts, nonz1.ecounts)
all.ecounts.z1 <- all.ecounts.z1[order(all.ecounts.z1[,1], decreasing = T),]
# Candidate 6979681 looks most promising, but it still has 11 apparent false positives. Looking at the clipped
# sequences shows that only two of these apparent false positives have the "right" clipped reads after the 
# breakpoint (AB0173_C and AB0178_C). But if you plot that sample, it turns out it is actually a true positive 
# that the coverage-based method missed. So it looks like this breakpoint is correct. This contrasts with the 
# other two candidates, where the clipped bases are the same between the dup samples and the false positives.

# So let's move forward with those two breakpoints
assumed.Dupz1.start.position <- 6968962
assumed.Dupz1.end.position <- 6979681
all.Dupz1.start.sequences <- character()
for (sample.name in sub('-', '_', Dupz1.samples)){
	these.start.sequences <- clipping.end.point.cyp.list[[sample.name]][clipping.end.point.cyp.list[[sample.name]]$Position == assumed.Dupz1.start.position, 'Clipped_sequence']
	all.Dupz1.start.sequences <- c(all.Dupz1.start.sequences, these.start.sequences)
}
longest.Dupz1.start.sequence <- all.Dupz1.start.sequences[which.max(nchar(all.Dupz1.start.sequences))]
all.Dupz1.end.sequences <- character()
for (sample.name in sub('-', '_', Dupz1.samples)){
	these.end.sequences <- clipping.start.point.cyp.list[[sample.name]][clipping.start.point.cyp.list[[sample.name]]$Position == assumed.Dupz1.end.position, 'Clipped_sequence']
	all.Dupz1.end.sequences <- c(all.Dupz1.end.sequences, these.end.sequences)
}
longest.Dupz1.end.sequence <- all.Dupz1.end.sequences[which.max(nchar(all.Dupz1.end.sequences))]
# The breakpoint looks like this (the CTCC is present at both ends and therefore could be either side of the breakpoint)
# TCGAATCAATCTTCGCA CTCC AGGTTAAACAACAGCGG
#  end of the dup ^      ^ start of the dup
#                 ^      ^ 
#  position 6979676      position 6968967

# What about CYP6M2?
# First, it looks like there are some cross-chrom reads that line up quite close to the end of the duplication.
# However, comparing these between samples shows that they are not consistent. Let's just look for shared 
# breakpoint reads instead. 
clipping.end.point.allm2dups <- do.call(rbind, clipping.end.point.list[sub('-', '_', encompass.duplication.detected.in.cyp6m2)])
clipping.end.point.allm2dups <- clipping.end.point.allm2dups[order(clipping.end.point.allm2dups$Position), ]
clipping.start.point.allm2dups <- do.call(rbind, clipping.start.point.list[sub('-', '_', encompass.duplication.detected.in.cyp6m2)])
clipping.start.point.allm2dups <- clipping.start.point.allm2dups[order(clipping.start.point.allm2dups$Position), ]
# Lets look at the potential start point. Because the region is so messy, we have little idea where the start 
# point could be, so we cast our net wide
clipping.end.point.allm2dups.narrow <- clipping.end.point.allm2dups[(clipping.end.point.allm2dups$Position > 6915000) & (clipping.end.point.allm2dups$Position < 6935000), ]
clipping.start.point.allm2dups.narrow <- clipping.start.point.allm2dups[(clipping.start.point.allm2dups$Position > 6915000) & (clipping.start.point.allm2dups$Position < 6935000), ]
plot(c(6915000, 6935000), c(0, 0), ylim = c(0,6), type = 'n')
add.breakpoints(clipping.end.point.allm2dups.narrow[, 'Position'], this.col = 'brown', yrange = c(0,6))
add.breakpoints(clipping.start.point.allm2dups.narrow[, 'Position'], this.col = 'pink', yrange = c(0,6))
abline(v = 6927000)
abline(v = cyp6m2.coords[2])
# That plot suggest some possible candidates: 6927828 for start point and 6930568 for end point
startpostable.m2 <- table(clipping.end.point.allm2dups.narrow$Position)
startpostable.m2 <- startpostable.m2[order(startpostable.m2, decreasing = T)]
endpostable.m2 <- table(clipping.start.point.allm2dups.narrow$Position)
endpostable.m2 <- endpostable.m2[order(endpostable.m2, decreasing = T)]
# Let's compare the frequency of these positions between individuals that do or don't have the dup (according
# to coverage).
# First the start point
m2 <- sub('-', '_', encompass.duplication.detected.in.cyp6m2)
nonm2 <- setdiff(names(counts.list), m2)
m2.scounts <- numeric()
nonm2.scounts <- numeric()
for (spos in names(startpostable.m2)){
	this.m2.count <- 0
	this.nonm2.count <- 0
	for (this.sample in m2)
		this.m2.count <- this.m2.count + any(clipping.end.point.cyp.list[[this.sample]]$Position == spos)
	for (this.sample in nonm2)
		this.nonm2.count <- this.nonm2.count + any(clipping.end.point.cyp.list[[this.sample]]$Position == spos)
	m2.scounts[as.character(spos)] <- this.m2.count
	nonm2.scounts[as.character(spos)] <- this.nonm2.count
}
all.scounts.m2 <- cbind(m2.scounts, nonm2.scounts)
all.scounts.m2 <- all.scounts.m2[order(all.scounts.m2[,1], decreasing = T),]
# 6927942 looks like it could be a candidate, though there are two false negatives. One of these false negatives
# is actually probably a true negative, since it looks like a different duplication, or possible not a duplication
# at all (AR0066_C).
# Now the end point. 
m2.ecounts <- numeric()
nonm2.ecounts <- numeric()
for (epos in names(endpostable.m2)){
	this.m2.count <- 0
	this.nonm2.count <- 0
	for (this.sample in m2)
		this.m2.count <- this.m2.count + any(clipping.start.point.list[[this.sample]]$Position == epos)
	for (this.sample in nonm2)
		this.nonm2.count <- this.nonm2.count + any(clipping.start.point.list[[this.sample]]$Position == epos)
	m2.ecounts[as.character(epos)] <- this.m2.count
	nonm2.ecounts[as.character(epos)] <- this.nonm2.count
}
all.ecounts.m2 <- cbind(m2.ecounts, nonm2.ecounts)
all.ecounts.m2 <- all.ecounts.m2[order(all.ecounts.m2[,1], decreasing = T),]
# Nothing looks great. The best is 6928838, with 3 false negatives and only 14 false positives. But the clipped
# bases are the same for the positives and negatives, so these are genuince false positives.
assumed.Dupm2.start.position <- 6927942
all.Dupm2.start.sequences <- character()
for (sample.name in sub('-', '_', Dupm2.samples)){
	these.start.sequences <- clipping.end.point.cyp.list[[sample.name]][clipping.end.point.cyp.list[[sample.name]]$Position == assumed.Dupm2.start.position, 'Clipped_sequence']
	all.Dupm2.start.sequences <- c(all.Dupm2.start.sequences, these.start.sequences)
}
longest.Dupm2.start.sequence <- all.Dupm2.start.sequences[which.max(nchar(all.Dupm2.start.sequences))]
# The breakpoint looks like this 
# TATATTATACCAAATTATTA   TACAAACCAAATTATACAAAA
#    end of the dup? ^   ^ start of the dup
#                    ^   ^ 
#     position ???????   position 6927943
# aligning the clipped bases before the known posiiton (6927942) to the genome mostly aligns around 30944010
# on 2L, (but only about 20 bases actually align) except for the longest string of bases, which aligns to 
# 13289520 on UNKN. 

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
	return.vector <- rep(2,ncol(read.based.cyp.duplications))
	names(return.vector) <- colnames(read.based.cyp.duplications)
	# Manually determine genomic regions for each of the duplication types. For Dupz1, the coverage increase
	# actually covers a smaller region than the breakpoints would suggest, so we give it a slightly smaller 
	# size
	Dupz1.region <- c(6970200, 6978000)
	Dupm2.region <- c(Dupm2.bp, 6930500)
	# Create regions as logical vectors for each of those
	pos.list <- list()
	pos.list[['Dupz1']] <- (counts.table$Position >= Dupz1.region[1]) & (counts.table$Position <= Dupz1.region[2] - coverage.window)
	pos.list[['Dupm2']] <- (counts.table$Position >= Dupm2.region[1]) & (counts.table$Position <= Dupm2.region[2] - coverage.window)
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
	sample.dups <- colnames(read.based.cyp.duplications)[read.based.cyp.duplications[sample.name, ] >= 1]
	sample.dups <- sample.dups[sample.dups != 'Dup0']
	remaining.dups <- sample.dups
	analysed.dups <- character()
	# If we have a Dup0, then there is something going on that we don't have a handle on so we bail
	Dup0 <- read.based.cyp.duplications[sample.name, 'Dup0']
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

# Let's have a look at the distribution of Dupz1
Dupz1.pop.table <- table(genotype.calls[,'Dupz1'], meta$population)
HWExact(c(as.numeric(Dupz1.pop.table[,'BFgam']), 0)) # We use "as.numeric" otherwise the function complains about the names

# Dupm2:
Dupm2.pop.table <- table(genotype.calls[,'Dupm2'], meta$population)
HWExact(c(as.numeric(Dupm2.pop.table[,'CIcol']), 0)) # We use "as.numeric" otherwise the function complains about the names

write.table(genotype.calls, 'genotype_calls.csv', sep = '\t', col.names = NA)

save.image('CYP6M2-Z1_analysis_shrunk_data.Rdata')



