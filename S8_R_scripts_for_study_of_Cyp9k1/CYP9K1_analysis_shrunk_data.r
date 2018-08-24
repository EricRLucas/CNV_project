library(stringr)
library(Biostrings)
library(HardyWeinberg)

load('CNV_analysis_phase2_for_CYP9K1.Rdata')

cyp9k1.coord <- c(15240572, 15240572 + 2292)
cyp9k1.region.plotting.indices <- c(min(which(counts.list[[1]]$Position >= cyp9k1.region.coord[1])), max(which(counts.list[[1]]$Position <= cyp9k1.region.coord[2])))

# Identify the samples with a detected duplication in cyp9k1
encompass.duplication.detected.in.cyp9k1 <- encompass.output.by.gene[['AGAP000818']]

# Load table containing geographical id
meta <- read.table('/home/eric/Liverpool/phase2.AR1/samples/samples.meta.txt', header = T, row.names = 1, sep = '\t', quote = "", comment.char = "")
rownames(meta) <- sub('-', '_', rownames(meta))
# Create booleans vector that tell us whether each individual in the samples table has a duplication in cyp9k1
hasencomp <- rownames(meta) %in% encompass.duplication.detected.in.cyp9k1
# add these data to the meta table
meta$hasencomp = hasencomp

# Now, within each population, we want to know the proportion of individuals carrying a duplication
propencomp <- tapply(meta$hasencomp, meta$population, mean)
pop.sizes <- tapply(meta$population, meta$population, length)
propdup.for.plot <- rbind(propencomp, 1 - propencomp)

x11(width = 15)
par(mfrow = c(1,1), mar = c(2,4,4,8), oma = c(3,0,0,0), xpd = F)
barplot(propdup.for.plot, ylab = 'Proportion', names.arg = rep('', ncol(propdup.for.plot)),  main = 'CYP9K1')
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
SSFA.folder <- paste('/home/eric/Liverpool/CNV_/SSFA/v3_', chrom, '/phase2/CYP9K1_region', sep = '')
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
FA.cyp9k1.list <- list()
for (i in 1:length(FA.list)){
	this.list <- FA.list[[i]]
	this.name <- names(FA.list)[i]
	# Use the following line if you want reads that start OR end in the region
	which.in.cyp9k1.FA <- ((this.list$Position >= cyp9k1.region.coord[1]) & (this.list$Position <= cyp9k1.region.coord[2])) | ((this.list$Mate.position >= cyp9k1.region.coord[1]) & (this.list$Mate.position <= cyp9k1.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the region
#	which.in.cyp9k1.FA <- ((this.list$Position >= cyp9k1.region.coord[1]) & (this.list$Mate.position <= cyp9k1.region.coord[2]))
	FA.cyp9k1.list[[this.name]] <- this.list[which.in.cyp9k1.FA, ]
}
SS.cyp9k1.list <- list()
for (i in 1:length(SS.list)){
	this.list <- SS.list[[i]]
	this.name <- names(SS.list)[i]
	# Use the following line if you want reads that start OR end in the region
	which.in.cyp9k1.SS <- ((this.list$Position >= cyp9k1.region.coord[1]) & (this.list$Position <= cyp9k1.region.coord[2])) | ((this.list$Mate.position >= cyp9k1.region.coord[1]) & (this.list$Mate.position <= cyp9k1.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the region
#	which.in.cyp9k1.SS <- ((this.list$Position >= cyp9k1.region.coord[1]) & (this.list$Mate.position <= cyp9k1.region.coord[2]))
	SS.cyp9k1.list[[this.name]] <- this.list[which.in.cyp9k1.SS, ]
}
FM.cyp9k1.list <- list()
for (i in 1:length(FM.list)){
	this.list <- FM.list[[i]]
	this.name <- names(FM.list)[i]
	# Use the following line if you want reads that start OR end in the region
	which.in.cyp9k1.FM <- ((this.list$Position >= cyp9k1.region.coord[1]) & (this.list$Position <= cyp9k1.region.coord[2])) | ((this.list$Mate.position >= cyp9k1.region.coord[1]) & (this.list$Mate.position <= cyp9k1.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the region
#	which.in.cyp9k1.FM <- ((this.list$Position >= cyp9k1.region.coord[1]) & (this.list$Mate.position <= cyp9k1.region.coord[2]))
	FM.cyp9k1.list[[this.name]] <- this.list[which.in.cyp9k1.FM, ]
}
crosschrom.cyp9k1.list <- list()
for (i in 1:length(crosschrom.list)){
	this.list <- crosschrom.list[[i]]
	this.name <- names(crosschrom.list)[i]
	which.in.cyp9k1.crosschrom <- (this.list$Position >= cyp9k1.region.coord[1]) & (this.list$Position <= cyp9k1.region.coord[2])
	crosschrom.cyp9k1.list[[this.name]] <- this.list[which.in.cyp9k1.crosschrom, ]
}

# Load up the results of the breakpoint detection 
cat('Loading CYP9K1 region breakpoints data\n')
bp.folder <- paste('/home/eric/Liverpool/CNV_/breakpoints/' , chrom, '/phase2/CYP9K1_region', sep = '')
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
clipping.start.point.cyp9k1.list <- list()
for (i in 1:length(clipping.start.point.list)){
	this.list <- clipping.start.point.list[[i]]
	this.name <- names(clipping.start.point.list)[i]
	which.in.cyp9k1.clipping.start.point <- (this.list$Position >= cyp9k1.region.coord[1]) & (this.list$Position <= cyp9k1.region.coord[2])
	this.initial.list <- this.list[which.in.cyp9k1.clipping.start.point, , drop = F]
	# The following for lines allow breakpoints to be kept only if their frequency exceeds a certain threshold. They 
	# are effectively pointless unless you set the threshold below above 0.
	reads.per.start.point <- tapply(this.initial.list$Position, as.factor(this.initial.list$Position), length)
	positions.to.keep <- this.initial.list$Position %in% names(reads.per.start.point)[reads.per.start.point > 0]
	clipping.start.point.cyp9k1.list[[this.name]] <- this.initial.list[positions.to.keep, ]
}
clipping.end.point.cyp9k1.list <- list()
for (i in 1:length(clipping.end.point.list)){
	this.list <- clipping.end.point.list[[i]]
	this.name <- names(clipping.end.point.list)[i]
	which.in.cyp9k1.clipping.end.point <- (this.list$Position >= cyp9k1.region.coord[1]) & (this.list$Position <= cyp9k1.region.coord[2])
	this.initial.list <- this.list[which.in.cyp9k1.clipping.end.point, , drop = F]
	reads.per.end.point <- tapply(this.initial.list$Position, as.factor(this.initial.list$Position), length)
	positions.to.keep <- this.initial.list$Position %in% names(reads.per.end.point)[reads.per.end.point > 0]
	clipping.end.point.cyp9k1.list[[this.name]] <- this.initial.list[positions.to.keep, ]
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
mapq0.proportions <- subset(mapq0.proportions, (Position >= cyp9k1.region.coord[1]) & (Position <= cyp9k1.region.coord[2]))
total.coverage <- mapq0.proportions$'Count mapq > 0' + mapq0.proportions$'Count mapq = 0'
mapq0.proportions$mapq0.prop <- mapq0.proportions$'Count mapq = 0' / total.coverage
# In some cases, there is no coverage, so we get infinite values. Let's deal with those
mapq0.proportions$mapq0.prop[total.coverage == 0] <- 0

# Create a list of colours that will be used for plotting the duplication types
duplication.colours <- c(rgb(0.6509804, 0.8078431, 0.8901961), rgb(0.1215686, 0.4705882, 0.7058824), rgb(0.6980392, 0.8745098, 0.5411765), rgb(0.2, 0.627451, 0.172549), rgb(0.9843137, 0.6039216, 0.6), rgb(0.8901961, 0.1019608, 0.1098039), rgb(0.9921569, 0.7490196, 0.4352941), rgb(1, 0.4980392, 0), rgb(0.7921569, 0.6980392, 0.8392157), rgb(0.4156863, 0.2392157, 0.6039216), rgb(1,0,1,0.7), rgb(0.6941176, 0.3490196, 0.1568627), rgb(0.5,0.5,0,0.8), rgb(0,0.5,0.5,0.8), rgb(0.5, 0, 0.5, 0.8), rgb(0.5, 0, 0))

plot.cyp9k1 <- function(this.sample, list.of.dups = NULL, FA = T, SS=T, FM=T, XC=T, BP=T, start.index = cyp9k1.region.plotting.indices[1], end.index = cyp9k1.region.plotting.indices[2]){
	this.sample.name <- regmatches(this.sample, regexpr('A.\\d\\d\\d\\d_C', this.sample))
	plot(counts.list[[this.sample]]$Position[start.index : end.index], counts.list[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		if (list.of.dups['Dup16'] >= 1){
			rect(15222810, -0.6, Dup16.end.bp, 0, col = duplication.colours[16], border = duplication.colours[16])
			text(Dup16.end.bp, -0.4, 'Dup16', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup15'] >= 1){
			rect(15233807, -0.6, Dup15.end.bp, 0, col = duplication.colours[15], border = duplication.colours[15])
			text(Dup15.end.bp, -0.4, 'Dup15', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup14'] >= 1){
			rect(15233807, -0.6, Dup14.end.bp, 0, col = duplication.colours[14], border = duplication.colours[14])
			text(Dup14.end.bp, -0.4, 'Dup14', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup13'] >= 1){
			rect(Dup13.bp[1], -0.6, Dup13.bp[2], 0, col = duplication.colours[13], border = duplication.colours[13])
			text(Dup13.bp[2], -0.4, 'Dup13', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup12'] >= 1){
			rect(Dup12.bp[1], -0.6, Dup12.bp[2], 0, col = duplication.colours[12], border = duplication.colours[12])
			text(Dup12.bp[2], -0.4, 'Dup12', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup11'] >= 1){
			rect(Dup11.bp[1], -0.6, Dup11.bp[2], 0, col = duplication.colours[11], border = duplication.colours[11])
			text(Dup11.bp[2], -0.4, 'Dup11', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup10'] >= 1){
			rect(Dup10.bp[1], -0.6, Dup10.bp[2], 0, col = duplication.colours[10], border = duplication.colours[10])
			text(Dup10.bp[2], -0.4, 'Dup10', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup9'] >= 1){
			rect(Dup9.bp[1], -0.6, Dup9.bp[2], 0, col = duplication.colours[9], border = duplication.colours[9])
			text(Dup9.bp[2], -0.4, 'Dup9', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup8'] >= 1){
			rect(Dup8.bp[1], -0.6, Dup8.bp[2], 0, col = duplication.colours[8], border = duplication.colours[8])
			text(Dup8.bp[2], -0.4, 'Dup8', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup7'] >= 1){
			rect(15238000, -0.6, 15246300,  0, col = duplication.colours[7], border = duplication.colours[7])
			text(15246300, -0.4, 'Dup7', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup6'] >= 1){
			rect(Dup6.bp[1], -0.6, Dup6.bp[2], 0, col = duplication.colours[6], border = duplication.colours[6])
			text(Dup6.bp[2], -0.4, 'Dup6', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup5'] >= 1){
			rect(Dup5.bp[1], -0.6, Dup5.bp[2], 0, col = duplication.colours[5], border = duplication.colours[5])
			text(Dup5.bp[2], -0.4, 'Dup5', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup4'] >= 1){
			rect(Dup4.bp[1], -0.6, Dup4.bp[2], 0, col = duplication.colours[4], border = duplication.colours[4])
			text(Dup4.bp[2], -0.4, 'Dup4', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup3'] >= 1){
			rect(15240464, -0.6, Dup3.end.bp, 0, col = duplication.colours[3], border = duplication.colours[3])
			text(Dup3.end.bp, -0.4, 'Dup3', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup2'] >= 1){
			rect(Dup2.bp[1], -0.6, Dup2.bp[2], 0, col = duplication.colours[2], border = duplication.colours[2])
			text(Dup2.bp[2], -0.4, 'Dup2', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup1'] >= 1){
			rect(Dup1.bp[1], -0.6, Dup1.bp[2], 0, col = duplication.colours[1], border = duplication.colours[1])
			text(Dup1.bp[2], -0.4, 'Dup1', cex = 1.2, adj = c(1,0))
		}
	}
	if (FA){
		if (nrow(FA.cyp9k1.list[[this.sample.name]]) > 0)
			add.faceaways(FA.cyp9k1.list[[this.sample.name]], this.col = 'blue', yrange = c(0,3))
	}
	if (SS){
		if (nrow(SS.cyp9k1.list[[this.sample.name]]) > 0)
			add.faceaways(SS.cyp9k1.list[[this.sample.name]], this.col = 'cyan', yrange = c(3,6))
	}
	if (FM){
		if (nrow(FM.cyp9k1.list[[this.sample.name]]) > 0)
			add.faceaways(FM.cyp9k1.list[[this.sample.name]], this.col = 'green', yrange = c(6,9))
	}
	if (XC){
		if (nrow(crosschrom.cyp9k1.list[[this.sample.name]])) 
			add.breakpoints(crosschrom.cyp9k1.list[[this.sample.name]]$Position, this.col = 'red', yrange = c(9,12))
	}
	if (BP){
		if (nrow(clipping.end.point.cyp9k1.list[[this.sample.name]]))
			add.breakpoints(clipping.end.point.cyp9k1.list[[this.sample.name]][, 'Position'], this.col = 'brown', yrange = c(0,12))
		if (nrow(clipping.start.point.cyp9k1.list[[this.sample.name]]))
			add.breakpoints(clipping.start.point.cyp9k1.list[[this.sample.name]][, 'Position'], this.col = 'pink', yrange = c(0,12))
	}
	these.cnv.states <- counts.list[[this.sample]]$CNV[start.index : end.index]
	lines(counts.list[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2)
	abline(v = cyp9k1.coord[1])
	abline(v = cyp9k1.coord[2])
	text(mean(cyp9k1.coord), 10, 'CYP9K1', srt=90, adj = 0)
}

plot.all.cyp9k1 <- function(list.of.all, list.of.dupl = encompass.duplication.detected.in.cyp9k1, FA=T, SS=T, FM=T, XC=T, BP=T, start.index = cyp9k1.region.plotting.indices[1], end.index = cyp9k1.region.plotting.indices[2]){
	x11()
	list.of.all <- sub('-', '_', list.of.all)
	list.of.dupl <- sub('-', '_', list.of.dupl)
	for (this.sample in list.of.all){
		if (this.sample %in% list.of.dupl)
			par(bg = 'beige')
		else
			par(bg = 'white')
		plot.cyp9k1(this.sample, NULL, FA, SS, FM, XC, BP, start.index, end.index)
		locator(1)
	}
}

plot.all.cyp9k1.with.read.dups <- function(list.of.samples, matrix.of.read.dups = read.based.cyp9k1.duplications, FA = T, SS = T, FM = T, XC=T, BP=T, start.index = cyp9k1.region.plotting.indices[1], end.index = cyp9k1.region.plotting.indices[2]){
	x11()
	x.midpoint <- mean(c(counts.list[[1]]$Position[start.index], counts.list[[1]]$Position[end.index]))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- sub('-', '_', list.of.samples[i])
		plot.cyp9k1(this.sample, matrix.of.read.dups[this.sample,], FA, SS, FM, XC, BP, start.index, end.index)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}


# Use the following definitions for the different duplications (numbers represent start range (1st row) and end range
# (2nd row) of the read pairs). These values were determined after visual inspection of the data. 
Dup1.fa = matrix(c(15242500, 15244500, 15242800, 15244800), 2, 2)
Dup1.bp = c(15242505,15244812) 
Dup1.bp.seq = c('GTTTG', 'CATAT')
#
Dup2.fa = matrix(c(15238300, 15240800, 15238600, 15241100), 2, 2)
Dup2.bp = c(15238400, 15241082)
Dup2.bp.seq = c('CCGGC',' CGGTA')
#
Dup3.fa = matrix(c(15240300, 15243450, 15240600, 15243750), 2, 2)
Dup3.end.bp = 15243860
Dup3.end.bp.seq = 'TGAAC'
#
Dup4.fa = matrix(c(15240600, 15244200, 15240900, 15244500), 2, 2)
Dup4.bp = c(15240608, 15244503)
Dup4.bp.seq = c('ATAAA', 'ACTGG')
#
Dup5.fa = matrix(c(15238800, 15243850, 15239100, 15244150), 2, 2)
Dup5.bp = c(15238911, 15244175)
Dup5.bp.seq = c('CACGT', 'AGTAA')
#
Dup6.fa = matrix(c(15236400, 15243250, 15236700, 15243550), 2, 2)
Dup6.bp = c(15236449, 15243646)
Dup6.bp.seq = c('TTTTT', 'GTTTT')
#
Dup7.ss = matrix(c(15245400, 15246900,15245700, 15247200), 2, 2)
Dup7.bp = c(15245768, 15247258)
Dup7.bp.seq <- c('TTTGT', 'TCTAA')
#
Dup8.fa = matrix(c(15239200, 15247250, 15239500, 15247550), 2, 2)
Dup8.bp = c(15239276, 15247645)
Dup8.bp.seq = c('AACAT', 'TTGCT')
#
Dup9.fa = matrix(c(15239100, 15248900, 15239400, 15249200), 2, 2)
Dup9.bp = c(15239184, 15249314)
Dup9.bp.seq = c('GCACA', 'AGTAC')
#
Dup10.fa = matrix(c(15234900, 15244750, 15235200, 15245050), 2, 2)
Dup10.bp = c(15234989, 15245128)
Dup10.bp.seq = c('GCACC', 'ATTCT')
#
Dup11.fa = matrix(c(15236900, 15246800, 15237200, 15247100), 2, 2)
Dup11.bp = c(15236922, 15247159)
Dup11.bp.seq = c('CATTA', 'TATCT')
#
Dup12.fa = matrix(c(15234400, 15244350, 15234700, 15244650), 2, 2)
Dup12.bp = c(15234434, 15244702)
Dup12.bp.seq = c('AACAG', 'TACTA')
#
Dup13.fa = matrix(c(15240100, 15250250, 15240400, 15250550), 2, 2)
Dup13.bp = c(15240067, 15250575)
Dup13.bp.seq = c('CCTAA', 'GTGTA')
#
# Dup 14 seems to have a different endpoint to Dup15, but the same insertion point way upstream
Dup14.fa = matrix(c(15244200, 9676400, 15244500, 9676700), 2, 2)
# Because Dup14 has the same start pos as dup15, we don't use the start breakpoint as a diagnostic
#Dup14.bp = c(15233807, 15244936)
#Dup14.bp.seq = c('GGGTT', 'CCCAA')
Dup14.end.bp = c(15244936)
Dup14.end.bp.seq = c('CCCAA')
#
Dup15.fa = matrix(c(15246250, 9676400, 15246550, 9676700), 2, 2)
# Because Dup14 has the same start pos as Dup15, we don't use the start breakpoint as a diagnostic
#Dup15.bp = c(15233807, 15246640)
#Dup15.bp.seq = c('GGGTT', 'CCCAA')
Dup15.end.bp = c(15246640)
Dup15.end.bp.seq = c('CCCAA')
#
Dup16.fa = matrix(c(15222700, 15244450, 15223000, 15244750), 2, 2)
Dup16.end.bp = 15244755
Dup16.end.bp.seq = 'AAGTA'

# We will store the duplication type information in a matrix
read.based.cyp9k1.duplications <- matrix(0, length(counts.list), 17, dimnames = list(names(counts.list), c('Dup0', 'Dup1', 'Dup2', 'Dup3', 'Dup4', 'Dup5', 'Dup6', 'Dup7', 'Dup8', 'Dup9', 'Dup10', 'Dup11', 'Dup12', 'Dup13', 'Dup14', 'Dup15', 'Dup16')))
# Now go through each individual and determine which duplications they have based on their diagnostic reads. 
for (this.sample in names(counts.list)){
	# get the FA and SS reads encompassed in the general region
	these.FA <- FA.cyp9k1.list[[this.sample]]
	these.SS <- SS.cyp9k1.list[[this.sample]]
	these.CSP <- clipping.start.point.cyp9k1.list[[this.sample]]
	these.CEP <- clipping.end.point.cyp9k1.list[[this.sample]]
	# Count the number of supporting reads
	num.Dup1.fa.reads <- sum((these.FA$Position > Dup1.fa[1,1]) & these.FA$Position < Dup1.fa[1,2]
	                       & (these.FA$Mate.position > Dup1.fa[2,1]) & (these.FA$Mate.position < Dup1.fa[2,2]))
	Dup1.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup1.bp[1]]
	num.Dup1.start.reads <- sum(substr(reverse(Dup1.CEP.seq), 1, 5) == Dup1.bp.seq[1])
	Dup1.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup1.bp[2]]
	num.Dup1.end.reads <- sum(substr(Dup1.CSP.seq, 1, 5) == Dup1.bp.seq[2])
	num.Dup1.reads <- num.Dup1.fa.reads + num.Dup1.start.reads + num.Dup1.end.reads
	#
	num.Dup2.fa.reads <- sum((these.FA$Position > Dup2.fa[1,1]) & these.FA$Position < Dup2.fa[1,2]
	                       & (these.FA$Mate.position > Dup2.fa[2,1]) & (these.FA$Mate.position < Dup2.fa[2,2]))
	Dup2.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup2.bp[1]]
	num.Dup2.start.reads <- sum(substr(reverse(Dup2.CEP.seq), 1, 5) == Dup2.bp.seq[1])
	Dup2.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup2.bp[2]]
	num.Dup2.end.reads <- sum(substr(Dup2.CSP.seq, 1, 5) == Dup2.bp.seq[2])
	num.Dup2.reads <- num.Dup2.fa.reads + num.Dup2.start.reads + num.Dup2.end.reads
	#
	num.Dup3.fa.reads <- sum((these.FA$Position > Dup3.fa[1,1]) & these.FA$Position < Dup3.fa[1,2]
	                       & (these.FA$Mate.position > Dup3.fa[2,1]) & (these.FA$Mate.position < Dup3.fa[2,2]))
	Dup3.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup3.end.bp]
	num.Dup3.end.reads <- sum(substr(Dup3.CSP.seq, 1, 5) == Dup3.end.bp.seq)
	num.Dup3.reads <- num.Dup3.fa.reads + num.Dup3.end.reads
	#
	num.Dup4.fa.reads <- sum((these.FA$Position > Dup4.fa[1,1]) & these.FA$Position < Dup4.fa[1,2]
	                       & (these.FA$Mate.position > Dup4.fa[2,1]) & (these.FA$Mate.position < Dup4.fa[2,2]))
	Dup4.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup4.bp[1]]
	num.Dup4.start.reads <- sum(substr(reverse(Dup4.CEP.seq), 1, 5) == Dup4.bp.seq[1])
	Dup4.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup4.bp[2]]
	num.Dup4.end.reads <- sum(substr(Dup4.CSP.seq, 1, 5) == Dup4.bp.seq[2])
	num.Dup4.reads <- num.Dup4.fa.reads + num.Dup4.start.reads + num.Dup4.end.reads
	#
	num.Dup5.fa.reads <- sum((these.FA$Position > Dup5.fa[1,1]) & these.FA$Position < Dup5.fa[1,2]
	                       & (these.FA$Mate.position > Dup5.fa[2,1]) & (these.FA$Mate.position < Dup5.fa[2,2]))
	Dup5.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup5.bp[1]]
	num.Dup5.start.reads <- sum(substr(reverse(Dup5.CEP.seq), 1, 5) == Dup5.bp.seq[1])
	Dup5.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup5.bp[2]]
	num.Dup5.end.reads <- sum(substr(Dup5.CSP.seq, 1, 5) == Dup5.bp.seq[2])
	num.Dup5.reads <- num.Dup5.fa.reads + num.Dup5.start.reads + num.Dup5.end.reads
	#
	num.Dup6.fa.reads <- sum((these.FA$Position > Dup6.fa[1,1]) & these.FA$Position < Dup6.fa[1,2]
	                       & (these.FA$Mate.position > Dup6.fa[2,1]) & (these.FA$Mate.position < Dup6.fa[2,2]))
	Dup6.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup6.bp[1]]
	num.Dup6.start.reads <- sum(substr(reverse(Dup6.CEP.seq), 1, 5) == Dup6.bp.seq[1])
	Dup6.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup6.bp[2]]
	num.Dup6.end.reads <- sum(substr(Dup6.CSP.seq, 1, 5) == Dup6.bp.seq[2])
	num.Dup6.reads <- num.Dup6.fa.reads + num.Dup6.start.reads + num.Dup6.end.reads
	#
	num.Dup7.ss.reads <- sum((these.SS$Position > Dup7.ss[1,1]) & these.SS$Position < Dup7.ss[1,2]
	                       & (these.SS$Mate.position > Dup7.ss[2,1]) & (these.SS$Mate.position < Dup7.ss[2,2]))
	Dup7.CSP.seq1 <- these.CSP$Clipped_sequence[these.CSP$Position == Dup7.bp[1]]
	num.Dup7.start.reads <- sum(substr(Dup7.CSP.seq1, 1, 5) == Dup7.bp.seq[1])
	Dup7.CSP.seq2 <- these.CSP$Clipped_sequence[these.CSP$Position == Dup7.bp[2]]
	num.Dup7.end.reads <- sum(substr(Dup7.CSP.seq2, 1, 5) == Dup7.bp.seq[2])
	num.Dup7.reads <- num.Dup7.ss.reads + num.Dup7.start.reads + num.Dup7.end.reads
	#
	num.Dup8.fa.reads <- sum((these.FA$Position > Dup8.fa[1,1]) & these.FA$Position < Dup8.fa[1,2]
	                      & (these.FA$Mate.position > Dup8.fa[2,1]) & (these.FA$Mate.position < Dup8.fa[2,2]))
	Dup8.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup8.bp[1]]
	num.Dup8.start.reads <- sum(substr(reverse(Dup8.CEP.seq), 1, 5) == Dup8.bp.seq[1])
	Dup8.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup8.bp[2]]
	num.Dup8.end.reads <- sum(substr(Dup8.CSP.seq, 1, 5) == Dup8.bp.seq[2])
	num.Dup8.reads <- num.Dup8.fa.reads + num.Dup8.start.reads + num.Dup8.end.reads
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
	num.Dup11.fa.reads <- sum((these.FA$Position > Dup11.fa[1,1]) & these.FA$Position < Dup11.fa[1,2]
	                       & (these.FA$Mate.position > Dup11.fa[2,1]) & (these.FA$Mate.position < Dup11.fa[2,2]))
	Dup11.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup11.bp[1]]
	num.Dup11.start.reads <- sum(substr(reverse(Dup11.CEP.seq), 1, 5) == Dup11.bp.seq[1])
	Dup11.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup11.bp[2]]
	num.Dup11.end.reads <- sum(substr(Dup11.CSP.seq, 1, 5) == Dup11.bp.seq[2])
	num.Dup11.reads <- num.Dup11.fa.reads + num.Dup11.start.reads + num.Dup11.end.reads
	#
	num.Dup12.fa.reads <- sum((these.FA$Position > Dup12.fa[1,1]) & these.FA$Position < Dup12.fa[1,2]
	                       & (these.FA$Mate.position > Dup12.fa[2,1]) & (these.FA$Mate.position < Dup12.fa[2,2]))
	Dup12.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup12.bp[1]]
	num.Dup12.start.reads <- sum(substr(reverse(Dup12.CEP.seq), 1, 5) == Dup12.bp.seq[1])
	Dup12.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup12.bp[2]]
	num.Dup12.end.reads <- sum(substr(Dup12.CSP.seq, 1, 5) == Dup12.bp.seq[2])
	num.Dup12.reads <- num.Dup12.fa.reads + num.Dup12.start.reads + num.Dup12.end.reads
	#
	num.Dup13.fa.reads <- sum((these.FA$Position > Dup13.fa[1,1]) & these.FA$Position < Dup13.fa[1,2]
	                       & (these.FA$Mate.position > Dup13.fa[2,1]) & (these.FA$Mate.position < Dup13.fa[2,2]))
	Dup13.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup13.bp[1]]
	num.Dup13.start.reads <- sum(substr(reverse(Dup13.CEP.seq), 1, 5) == Dup13.bp.seq[1])
	Dup13.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup13.bp[2]]
	num.Dup13.end.reads <- sum(substr(Dup13.CSP.seq, 1, 5) == Dup13.bp.seq[2])
	num.Dup13.reads <- num.Dup13.fa.reads + num.Dup13.start.reads + num.Dup13.end.reads
	#
	num.Dup14.fa.reads <- sum((these.FA$Position > Dup14.fa[1,1]) & these.FA$Position < Dup14.fa[1,2]
	                       & (these.FA$Mate.position > Dup14.fa[2,1]) & (these.FA$Mate.position < Dup14.fa[2,2]))
	Dup14.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup14.end.bp]
	num.Dup14.end.reads <- sum(substr(Dup14.CSP.seq, 1, 5) == Dup14.end.bp.seq)
	num.Dup14.reads <- num.Dup14.fa.reads + num.Dup14.end.reads
	#
	num.Dup15.fa.reads <- sum((these.FA$Position > Dup15.fa[1,1]) & these.FA$Position < Dup15.fa[1,2]
	                       & (these.FA$Mate.position > Dup15.fa[2,1]) & (these.FA$Mate.position < Dup15.fa[2,2]))
	Dup15.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup15.end.bp]
	num.Dup15.end.reads <- sum(substr(Dup15.CSP.seq, 1, 5) == Dup15.end.bp.seq)
	num.Dup15.reads <- num.Dup15.fa.reads + num.Dup15.end.reads
	#
	num.Dup16.fa.reads <- sum((these.FA$Position > Dup16.fa[1,1]) & these.FA$Position < Dup16.fa[1,2]
	                       & (these.FA$Mate.position > Dup16.fa[2,1]) & (these.FA$Mate.position < Dup16.fa[2,2]))
	Dup16.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup16.end.bp]
	num.Dup16.end.reads <- sum(substr(Dup16.CSP.seq, 1, 5) == Dup16.end.bp.seq)
	num.Dup16.reads <- num.Dup16.fa.reads + num.Dup16.end.reads
	# 
	# For each of these, if there are at least 2 supporting reads, score it as a duplication. If there is only one
	# supporting read, score it tentatively
	read.based.cyp9k1.duplications[this.sample, 2:17] <- as.logical(c(num.Dup1.reads, num.Dup2.reads, num.Dup3.reads, num.Dup4.reads, num.Dup5.reads, num.Dup6.reads, num.Dup7.reads, num.Dup8.reads, num.Dup9.reads, num.Dup10.reads, num.Dup11.reads, num.Dup12.reads, num.Dup13.reads, num.Dup14.reads, num.Dup15.reads, num.Dup16.reads))
	if (num.Dup1.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup1'] <- 0.5
	if (num.Dup2.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup2'] <- 0.5
	if (num.Dup3.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup3'] <- 0.5
	if (num.Dup4.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup4'] <- 0.5
	if (num.Dup5.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup5'] <- 0.5
	if (num.Dup6.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup6'] <- 0.5
	if (num.Dup7.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup7'] <- 0.5
	if (num.Dup8.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup8'] <- 0.5
	if (num.Dup9.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup9'] <- 0.5
	if (num.Dup10.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup10'] <- 0.5
	if (num.Dup11.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup11'] <- 0.5
	if (num.Dup12.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup12'] <- 0.5
	if (num.Dup13.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup13'] <- 0.5
	if (num.Dup14.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup14'] <- 0.5
	if (num.Dup15.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup15'] <- 0.5
	if (num.Dup16.reads == 1)
		read.based.cyp9k1.duplications[this.sample, 'Dup16'] <- 0.5
	# If the sample has a duplication, but is none of the specific ones, we record it as Dup0
	if ((sum(read.based.cyp9k1.duplications[this.sample, ] >= 1) == 0) & (this.sample %in% encompass.duplication.detected.in.cyp9k1))
		read.based.cyp9k1.duplications[this.sample, 'Dup0'] <- 1
}

# Compared the duplications detected by coverage or reads
cov.based.duplications.names <- encompass.duplication.detected.in.cyp9k1
read.based.duplications.names <- rownames(read.based.cyp9k1.duplications)[as.logical(apply(read.based.cyp9k1.duplications[,2:17] >= 1, 1, sum))]

cov.based.negatives <- setdiff(read.based.duplications.names, cov.based.duplications.names)
read.based.negatives <- setdiff(cov.based.duplications.names, read.based.duplications.names)

# Output a table of each duplication type and the country they are found in
Dup0.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup0'] >= 1], 'population', drop = F]
Dup1.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup1'] >= 1], 'population', drop = F]
Dup2.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup2'] >= 1], 'population', drop = F]
Dup3.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup3'] >= 1], 'population', drop = F]
Dup4.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup4'] >= 1], 'population', drop = F]
Dup5.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup5'] >= 1], 'population', drop = F]
Dup6.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup6'] >= 1], 'population', drop = F]
Dup7.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup7'] >= 1], 'population', drop = F]
Dup8.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup8'] >= 1], 'population', drop = F]
Dup9.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup9'] >= 1], 'population', drop = F]
Dup10.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup10'] >= 1], 'population', drop = F]
Dup11.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup11'] >= 1], 'population', drop = F]
Dup12.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup12'] >= 1], 'population', drop = F]
Dup13.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup13'] >= 1], 'population', drop = F]
Dup14.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup14'] >= 1], 'population', drop = F]
Dup15.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup15'] >= 1], 'population', drop = F]
Dup16.meta <- meta[rownames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[,'Dup16'] >= 1], 'population', drop = F]

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
write(rownames(Dup12.meta), file = 'Dup12_samples.txt', ncolumns = 1)
write(rownames(Dup13.meta), file = 'Dup13_samples.txt', ncolumns = 1)
write(rownames(Dup14.meta), file = 'Dup14_samples.txt', ncolumns = 1)
write(rownames(Dup15.meta), file = 'Dup15_samples.txt', ncolumns = 1)
write(rownames(Dup16.meta), file = 'Dup16_samples.txt', ncolumns = 1)

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
Dup12.samples <- rownames(Dup12.meta)
Dup13.samples <- rownames(Dup13.meta)
Dup14.samples <- rownames(Dup14.meta)
Dup15.samples <- rownames(Dup15.meta)
Dup16.samples <- rownames(Dup16.meta)


######################################
# PRECISELY IDENTIFYING BREAK POINTS #
######################################

# Load the genome
genome <- readDNAStringSet('/home/eric/Liverpool/AR3/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa')

# Dup1:
assumed.Dup1.start.position <- 15242505
assumed.Dup1.end.position <- 15244812
all.Dup1.start.sequences <- character()
for (sample.name in sub('-', '_', Dup1.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup1.start.position, 'Clipped_sequence']
	all.Dup1.start.sequences <- c(all.Dup1.start.sequences, these.start.sequences)
}
longest.Dup1.start.sequence <- all.Dup1.start.sequences[which.max(nchar(all.Dup1.start.sequences))]
# Do the same at the end of the sequence
all.Dup1.end.sequences <- character()
for (sample.name in sub('-', '_', Dup1.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup1.end.position, 'Clipped_sequence']
	all.Dup1.end.sequences <- c(all.Dup1.end.sequences, these.end.sequences)
}
longest.Dup1.end.sequence <- all.Dup1.end.sequences[which.max(nchar(all.Dup1.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup1 <- substr(genome[6], assumed.Dup1.end.position-100, assumed.Dup1.end.position-1)
right.Dup1 <- substr(genome[6], assumed.Dup1.start.position+1, assumed.Dup1.start.position+100)
# Dup1 doesn't cover the whole of CYP9K1. The start point is near the end of the last exon. 
# The breakpoint looks like this 
# CATCGCTTTGCCGTTTG   CATATCGCTCAGCACATC
#  end of the dup ^   ^ start of the dup
#                 ^   ^
# position 15244811   position 15242506

# Dup2:
assumed.Dup2.start.position <- 15238400
assumed.Dup2.end.position <- 15241082
all.Dup2.start.sequences <- character()
for (sample.name in sub('-', '_', Dup2.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup2.start.position, 'Clipped_sequence']
	all.Dup2.start.sequences <- c(all.Dup2.start.sequences, these.start.sequences)
}
longest.Dup2.start.sequence <- all.Dup2.start.sequences[which.max(nchar(all.Dup2.start.sequences))]
# Do the same at the end of the sequence
all.Dup2.end.sequences <- character()
for (sample.name in sub('-', '_', Dup2.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup2.end.position, 'Clipped_sequence']
	all.Dup2.end.sequences <- c(all.Dup2.end.sequences, these.end.sequences)
}
longest.Dup2.end.sequence <- all.Dup2.end.sequences[which.max(nchar(all.Dup2.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup2 <- substr(genome[6], assumed.Dup2.end.position-100, assumed.Dup2.end.position-1)
right.Dup2 <- substr(genome[6], assumed.Dup2.start.position+1, assumed.Dup2.start.position+100)
# The clipped sequences from 15238400 align to 15241081 after an insertion that is mostly a duplicate of the sequence 
# before the breakpoint. The duplication ends about half-way through the first exon of CYP9K1. 
# The breakpoint looks like this
# GTCGCGCATGATGGCCGCGGCCTGTATCCAC CGGTAGCATGATGGCCGCGGCC GGTAGCCGCTTACTGTTGGT
#                end of the dup ^      inserted seq      ^ start of the dup
#                               ^                        ^
#               position 15241081                        position 15238401

# Dup3:
assumed.Dup3.end.position <- 15243860
all.Dup3.end.sequences <- character()
for (sample.name in sub('-', '_', Dup3.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup3.end.position, 'Clipped_sequence']
	all.Dup3.end.sequences <- c(all.Dup3.end.sequences, these.end.sequences)
}
longest.Dup3.end.sequence <- all.Dup3.end.sequences[which.max(nchar(all.Dup3.end.sequences))]
# The clipped sequences align to 15240464. This makes sense with the FA reads. There aren't any obvious reads at the end
# breakpoint in the one individual with this duplication, so can't get confirmation that this is correct. If the breakpoints
# have been correctly identified, then they are just before and just after CYP9K1. 
# The breakpoint looks like this 
# AGGAAACTAATATATGTATTTA      T     GAACTAAATCACCGGTTAGG
#       end of the dup ^  insertion ^ start of the dup
#                      ^            ^
#      position 15243859            position 15240464

# Dup4:
assumed.Dup4.start.position <- 15240608
assumed.Dup4.end.position <- 15244503
all.Dup4.start.sequences <- character()
for (sample.name in sub('-', '_', Dup4.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup4.start.position, 'Clipped_sequence']
	all.Dup4.start.sequences <- c(all.Dup4.start.sequences, these.start.sequences)
}
longest.Dup4.start.sequence <- all.Dup4.start.sequences[which.max(nchar(all.Dup4.start.sequences))]
# Do the same at the end of the sequence
all.Dup4.end.sequences <- character()
for (sample.name in sub('-', '_', Dup4.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup4.end.position, 'Clipped_sequence']
	all.Dup4.end.sequences <- c(all.Dup4.end.sequences, these.end.sequences)
}
longest.Dup4.end.sequence <- all.Dup4.end.sequences[which.max(nchar(all.Dup4.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup4 <- substr(genome[6], assumed.Dup4.end.position-100, assumed.Dup4.end.position-1)
right.Dup4 <- substr(genome[6], assumed.Dup4.start.position+1, assumed.Dup4.start.position+100)
# The breakpoint looks like this (the A is inserted between the two ends of the breakpoint)
# The breakpoint occurs in the first intron (before the first exon). 
# GTTATTAAGAATGTGAAAT        A       CTGGACACTAAACGCTCAA
#    end of the dup ^  inserted seq  ^ start of the dup
#                   ^                ^ 
#   position 15244502                position 15240609

# Dup5:
assumed.Dup5.start.position <- 15238911
assumed.Dup5.end.position <- 15244175
all.Dup5.start.sequences <- character()
for (sample.name in sub('-', '_', Dup5.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup5.start.position, 'Clipped_sequence']
	all.Dup5.start.sequences <- c(all.Dup5.start.sequences, these.start.sequences)
}
longest.Dup5.start.sequence <- all.Dup5.start.sequences[which.max(nchar(all.Dup5.start.sequences))]
# Do the same at the end of the sequence
all.Dup5.end.sequences <- character()
for (sample.name in sub('-', '_', Dup5.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup5.end.position, 'Clipped_sequence']
	all.Dup5.end.sequences <- c(all.Dup5.end.sequences, these.end.sequences)
}
longest.Dup5.end.sequence <- all.Dup5.end.sequences[which.max(nchar(all.Dup5.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup5 <- substr(genome[6], assumed.Dup5.end.position-100, assumed.Dup5.end.position-1)
right.Dup5 <- substr(genome[6], assumed.Dup5.start.position+1, assumed.Dup5.start.position+100)
# The breakpoint looks like this (the TA is present at both ends and therefore could be either side of the breakpoint)
# TTGGACCATCATTGCAC TA AGTAAAAATAATTGTA
#  end of the dup ^    ^ start of the dup
#                 ^    ^ 
# position 15244172    position 15238914

# Dup6:
assumed.Dup6.start.position <- 15236449
assumed.Dup6.end.position <- 15243646
all.Dup6.start.sequences <- character()
for (sample.name in sub('-', '_', Dup6.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup6.start.position, 'Clipped_sequence']
	all.Dup6.start.sequences <- c(all.Dup6.start.sequences, these.start.sequences)
}
longest.Dup6.start.sequence <- all.Dup6.start.sequences[which.max(nchar(all.Dup6.start.sequences))]
# Do the same at the end of the sequence
all.Dup6.end.sequences <- character()
for (sample.name in sub('-', '_', Dup6.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup6.end.position, 'Clipped_sequence']
	all.Dup6.end.sequences <- c(all.Dup6.end.sequences, these.end.sequences)
}
longest.Dup6.end.sequence <- all.Dup6.end.sequences[which.max(nchar(all.Dup6.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup6 <- substr(genome[6], assumed.Dup6.end.position-100, assumed.Dup6.end.position-1)
right.Dup6 <- substr(genome[6], assumed.Dup6.start.position+1, assumed.Dup6.start.position+100)
# The breakpoint looks like this (it has an unknown sequence inserted inside it)
# GATTGACGCAGTGGTTATC GTTTTTGCGGTTTTTT GTAGTACGTAGAGATGG
#    end of the dup ^  inserted seq    ^ start of the dup
#                   ^                  ^
#   position 15243646                  position 15236450
# I initially tried a different start point: 15236370. This looked good because it was close to the FA reads position,
# but this is actually an insertion that forces the reads to be clipped. It looks like this:
# TATAGACC CAGTTGCACTAT ATACAGTGCAAAGCCCGGTTAT
#        ^  insertion   ^
# 15236370              15236371

# Dup7 tends to be associated with an SS duplication to its right (which I'm not counting because it 
# doesn't cover CYP9K1). This SS duplication is defined by SS reads ~15245600 to ~15247150 and BP reads clipped at 
# 15245768 and 15247258 (both clipped at the end of the read, since the duplication is inverted). The clipped sequence 
# aligns where expected at the other breakpoint. There is also a deletion in some samples (again, doesn't cover CYP9K1
# so I don't consider it further) supported by FM reads ~15236150 to ~15237800 and BP reads clipped at 15236308 and
# 15237682. Again, the clipped sequences align where expected at the other breakpoint. 
# The duplication itself is harder to pin down. It may be associated with some BP reads that are clipped at 15245285. 
# The clipped sequence for this BP aligns equally at 9672726 or 9680423 (in opposite orientation to 9672726). I Don't 
# even know yet if this is a pragmatic solution, but it certainly doesn't get me any closer to understanding the 
# duplication (there are no FA or SS reads). I've tried implementing that diagnostic and I don't really believe it. 
# there are a lot of samples that it doesn't catch. 
# The Dup7 duplication looks like it could be a triplication where the third copy is inverted and has a small deletion.
# This would explain the SS reads (marking the deletion in the inverted copy) and why the coverage  in the region 
# marked by the SS reads tends to be 0.5, while the coverage in the rest of the region is 1. There is also a region
# to the right of the 1 region that seems to have coverage 0.5, but I'm not sure. 
# We used the SS reads and their associated breakpoints to diagnose Dup7. It leaves us not understanding the 
# breakpoint. I don't think it is just another manifestation of Dup15 because several Dup7 samples don't have Dup15.
assumed.Dup7.ss.start.position <- 15245768
assumed.Dup7.ss.end.position <- 15247258
all.Dup7.start.sequences <- character()
for (sample.name in sub('-', '_', Dup7.samples)){
	these.start.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup7.ss.start.position, 'Clipped_sequence']
	all.Dup7.start.sequences <- c(all.Dup7.start.sequences, these.start.sequences)
}
longest.Dup7.start.sequence <- all.Dup7.start.sequences[which.max(nchar(all.Dup7.start.sequences))]
# Do the same at the end of the sequence
all.Dup7.end.sequences <- character()
for (sample.name in sub('-', '_', Dup7.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup7.ss.end.position, 'Clipped_sequence']
	all.Dup7.end.sequences <- c(all.Dup7.end.sequences, these.end.sequences)
}
longest.Dup7.end.sequence <- all.Dup7.end.sequences[which.max(nchar(all.Dup7.end.sequences))]

# Dup8:
assumed.Dup8.start.position <- 15239276
assumed.Dup8.end.position <- 15247645
all.Dup8.start.sequences <- character()
for (sample.name in sub('-', '_', Dup8.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup8.start.position, 'Clipped_sequence']
	all.Dup8.start.sequences <- c(all.Dup8.start.sequences, these.start.sequences)
}
longest.Dup8.start.sequence <- all.Dup8.start.sequences[which.max(nchar(all.Dup8.start.sequences))]
# Do the same at the end of the sequence
all.Dup8.end.sequences <- character()
for (sample.name in sub('-', '_', Dup8.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup8.end.position, 'Clipped_sequence']
	all.Dup8.end.sequences <- c(all.Dup8.end.sequences, these.end.sequences)
}
longest.Dup8.end.sequence <- all.Dup8.end.sequences[which.max(nchar(all.Dup8.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup8 <- substr(genome[6], assumed.Dup8.end.position-100, assumed.Dup8.end.position-1)
right.Dup8 <- substr(genome[6], assumed.Dup8.start.position+1, assumed.Dup8.start.position+100)
# The breakpoint looks like this (the CGTA is present at both ends and therefore could be either side of the breakpoint)
# GCTGTTGCAATGCTACAA CGTA TTGCTTAAGAAAAGTA
#   end of the dup ^      ^ start of the dup
#                  ^      ^
#  position 15247640      position 15239281
# This duplication is also associated with a big deletion inside it (most of the second half of the duplication is 
# deleted), but this doesn't cover CYP9K1, it is supported by FM reads and BP reads at 15244041 and 15247267 (these 
# are in the right position relative to the FM reads, and the clipped sequences align in the right place). 

# Dup9:
assumed.Dup9.start.position <- 15239184
assumed.Dup9.end.position <- 15249314
all.Dup9.start.sequences <- character()
for (sample.name in sub('-', '_', Dup9.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup9.start.position, 'Clipped_sequence']
	all.Dup9.start.sequences <- c(all.Dup9.start.sequences, these.start.sequences)
}
longest.Dup9.start.sequence <- all.Dup9.start.sequences[which.max(nchar(all.Dup9.start.sequences))]
# Do the same at the end of the sequence
all.Dup9.end.sequences <- character()
for (sample.name in sub('-', '_', Dup9.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup9.end.position, 'Clipped_sequence']
	all.Dup9.end.sequences <- c(all.Dup9.end.sequences, these.end.sequences)
}
longest.Dup9.end.sequence <- all.Dup9.end.sequences[which.max(nchar(all.Dup9.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup9 <- substr(genome[6], assumed.Dup9.end.position-100, assumed.Dup9.end.position-1)
right.Dup9 <- substr(genome[6], assumed.Dup9.start.position+1, assumed.Dup9.start.position+100)
# The breakpoint looks like this. The GT is present at both ends and therefore could be either side of the breakpoint.
# The first 9 bases of at the start (GAACAGTGA) are reasonably close to the bases just after the end breakpoint (TAATAGAGA)
# which means that they actually align, hence the fact that the assumed.Dup9.end.position (where the reads are clipped) 
# is larger than the actual start point.
# GGCCTCCTAGATACACG GT GAACAGTGAAGTACTATTTTTTTCG
#  end of the dup ^    ^ start of the dup
#                 ^    ^ 
# position 15249302    position 15239187

# Dup10:
assumed.Dup10.start.position <- 15234989
assumed.Dup10.end.position <- 15245128
all.Dup10.start.sequences <- character()
for (sample.name in sub('-', '_', Dup10.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup10.start.position, 'Clipped_sequence']
	all.Dup10.start.sequences <- c(all.Dup10.start.sequences, these.start.sequences)
}
longest.Dup10.start.sequence <- all.Dup10.start.sequences[which.max(nchar(all.Dup10.start.sequences))]
# Do the same at the end of the sequence
all.Dup10.end.sequences <- character()
for (sample.name in sub('-', '_', Dup10.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup10.end.position, 'Clipped_sequence']
	all.Dup10.end.sequences <- c(all.Dup10.end.sequences, these.end.sequences)
}
longest.Dup10.end.sequence <- all.Dup10.end.sequences[which.max(nchar(all.Dup10.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup10 <- substr(genome[6], assumed.Dup10.end.position-100, assumed.Dup10.end.position-1)
right.Dup10 <- substr(genome[6], assumed.Dup10.start.position+1, assumed.Dup10.start.position+100)
# The breakpoint looks like this (the C is present at both ends and therefore could be either side of the breakpoint)
# ACACAGGCTCGACACATCCGCCACG C ATTCTGAAAAGCCTTATCAGTT
#          end of the dup ^   ^ start of the dup
#                         ^   ^
#         position 15245126   position 15234991

# Dup11:
assumed.Dup11.start.position <- 15236922
assumed.Dup11.end.position <- 15247159
all.Dup11.start.sequences <- character()
for (sample.name in sub('-', '_', Dup11.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup11.start.position, 'Clipped_sequence']
	all.Dup11.start.sequences <- c(all.Dup11.start.sequences, these.start.sequences)
}
longest.Dup11.start.sequence <- all.Dup11.start.sequences[which.max(nchar(all.Dup11.start.sequences))]
# Do the same at the end of the sequence
all.Dup11.end.sequences <- character()
for (sample.name in sub('-', '_', Dup11.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup11.end.position, 'Clipped_sequence']
	all.Dup11.end.sequences <- c(all.Dup11.end.sequences, these.end.sequences)
}
longest.Dup11.end.sequence <- all.Dup11.end.sequences[which.max(nchar(all.Dup11.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup11 <- substr(genome[6], assumed.Dup11.end.position-100, assumed.Dup11.end.position-1)
right.Dup11 <- substr(genome[6], assumed.Dup11.start.position+1, assumed.Dup11.start.position+100)
# A bunch of bases are inserted in the breakpoint. It looks like it's a repeat of the sequence just before the 
# breakpoint.
# AATAGTTATGTTTCTAGTTTTATTACATAT TATCTAGTTTTATTAC ATACGTGCTACCGAGTTTTAACCGT
#               end of the dup ^    unkown seq    ^ start of the dup
#                              ^                  ^
#              position 15247158                  position 15236923

# Dup12:
assumed.Dup12.start.position <- 15234434
assumed.Dup12.end.position <- 15244702
all.Dup12.start.sequences <- character()
for (sample.name in sub('-', '_', Dup12.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup12.start.position, 'Clipped_sequence']
	all.Dup12.start.sequences <- c(all.Dup12.start.sequences, these.start.sequences)
}
longest.Dup12.start.sequence <- all.Dup12.start.sequences[which.max(nchar(all.Dup12.start.sequences))]
# Do the same at the end of the sequence
all.Dup12.end.sequences <- character()
for (sample.name in sub('-', '_', Dup12.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup12.end.position, 'Clipped_sequence']
	all.Dup12.end.sequences <- c(all.Dup12.end.sequences, these.end.sequences)
}
longest.Dup12.end.sequence <- all.Dup12.end.sequences[which.max(nchar(all.Dup12.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup12 <- substr(genome[6], assumed.Dup12.end.position-100, assumed.Dup12.end.position-1)
right.Dup12 <- substr(genome[6], assumed.Dup12.start.position+1, assumed.Dup12.start.position+100)
# The breakpoint looks like this (the G is present at both ends and therefore could be either side of the breakpoint)
# TGGGAACTGCTTGTGACAA G TACTATATAGGTACT
#    end of the dup ^   ^ start of the dup
#                   ^   ^
#   position 15244700   position 15234436

# Dup13:
assumed.Dup13.start.position <- 15240067
assumed.Dup13.end.position <- 15250575
all.Dup13.start.sequences <- character()
for (sample.name in sub('-', '_', Dup13.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup13.start.position, 'Clipped_sequence']
	all.Dup13.start.sequences <- c(all.Dup13.start.sequences, these.start.sequences)
}
longest.Dup13.start.sequence <- all.Dup13.start.sequences[which.max(nchar(all.Dup13.start.sequences))]
# Do the same at the end of the sequence
all.Dup13.end.sequences <- character()
for (sample.name in sub('-', '_', Dup13.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup13.end.position, 'Clipped_sequence']
	all.Dup13.end.sequences <- c(all.Dup13.end.sequences, these.end.sequences)
}
longest.Dup13.end.sequence <- all.Dup13.end.sequences[which.max(nchar(all.Dup13.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
# At the left of the start point, we expect the reads to be the left of the end point. So let's find the 100bp to
# the left of position 28598850
left.Dup13 <- substr(genome[6], assumed.Dup13.end.position-100, assumed.Dup13.end.position-1)
right.Dup13 <- substr(genome[6], assumed.Dup13.start.position+1, assumed.Dup13.start.position+100)
# The breakpoint looks like this (the GA is present at both ends and therefore could be either side of the breakpoint)
# TTAGACATGGTAAAATCC GA GTGTAAGGCAATTTTAA
#   end of the dup ^    ^ start of the dup
#                  ^    ^
#  position 15250572    position 15240070

# Dup14:
assumed.Dup14.start.position <- 15233807 # This is exactly the same as the Dup15 start point. A coincidence?
assumed.Dup14.end.position <- 15244936
all.Dup14.start.sequences <- character()
for (sample.name in sub('-', '_', Dup14.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup14.start.position, 'Clipped_sequence']
	all.Dup14.start.sequences <- c(all.Dup14.start.sequences, these.start.sequences)
}
longest.Dup14.start.sequence <- all.Dup14.start.sequences[which.max(nchar(all.Dup14.start.sequences))]
# Do the same at the end of the sequence
all.Dup14.end.sequences <- character()
for (sample.name in sub('-', '_', Dup14.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup14.end.position, 'Clipped_sequence']
	all.Dup14.end.sequences <- c(all.Dup14.end.sequences, these.end.sequences)
}
longest.Dup14.end.sequence <- all.Dup14.end.sequences[which.max(nchar(all.Dup14.end.sequences))]
# The trailing regions are reverse complements of each-other:
clip1 <- longest.Dup14.end.sequence
clip2 <- reverseComplement(DNAString(longest.Dup14.start.sequence))
# This duplication has a different end point to Dup15, but its FA reads map to the same position ~5 Mbp upstream.
# Sure enough, the BP reads that map to the end of the duplication have their clipped sequences mapped to exactly 
# the same place as Dup15 (either before 9676356 or 9676537, since these two options are identical). 
# The FA reads stretch from the end of the dup to somewhere ~5 Mbp upstream.
# Oddly, the clipped characters either side of the breakpoints are reverse complements of each other, and map to the
# same 5 Mbp upstream region that the FA reads map to. There are two such mapping sites, close to each other and 
# perfect reverse-complements of each-other (9676277-9676356 and 9676537-9676616). 

# Dup15:
assumed.Dup15.start.position <- 15233807
assumed.Dup15.end.position <- 15246640
all.Dup15.start.sequences <- character()
for (sample.name in sub('-', '_', Dup15.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup15.start.position, 'Clipped_sequence']
	all.Dup15.start.sequences <- c(all.Dup15.start.sequences, these.start.sequences)
}
longest.Dup15.start.sequence <- all.Dup15.start.sequences[which.max(nchar(all.Dup15.start.sequences))]
# Do the same at the end of the sequence
all.Dup15.end.sequences <- character()
for (sample.name in sub('-', '_', Dup15.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup15.end.position, 'Clipped_sequence']
	all.Dup15.end.sequences <- c(all.Dup15.end.sequences, these.end.sequences)
}
longest.Dup15.end.sequence <- all.Dup15.end.sequences[which.max(nchar(all.Dup15.end.sequences))]
# The trailing regions are reverse complements of each-other:
clip1 <- longest.Dup15.end.sequence
clip2 <- reverseComplement(DNAString(longest.Dup15.start.sequence))
# The FA reads stretch from the end of the dup to somewhere ~5 Mbp upstream.
# Oddly, the clipped characters either side of the breakpoints are reverse complements of each other, and map to the
# same 5 Mbp upstream region that the FA reads map to. There are two such mapping sites, close to each other and 
# perfect reverse-complements of each-other (9676277-9676356 and 9676537-9676616). 

# Dup16:
assumed.Dup16.start.position <- 15222981
assumed.Dup16.end.position <- 15244755
all.Dup16.start.sequences <- character()
for (sample.name in sub('-', '_', Dup16.samples)){
	these.start.sequences <- clipping.end.point.cyp9k1.list[[sample.name]][clipping.end.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup16.start.position, 'Clipped_sequence']
	all.Dup16.start.sequences <- c(all.Dup16.start.sequences, these.start.sequences)
}
longest.Dup16.start.sequence <- all.Dup16.start.sequences[which.max(nchar(all.Dup16.start.sequences))]
# Do the same at the end of the sequence
all.Dup16.end.sequences <- character()
for (sample.name in sub('-', '_', Dup16.samples)){
	these.end.sequences <- clipping.start.point.cyp9k1.list[[sample.name]][clipping.start.point.cyp9k1.list[[sample.name]]$Position == assumed.Dup16.end.position, 'Clipped_sequence']
	all.Dup16.end.sequences <- c(all.Dup16.end.sequences, these.end.sequences)
}
longest.Dup16.end.sequence <- all.Dup16.end.sequences[which.max(nchar(all.Dup16.end.sequences))]
# Get the 100bp expected on either side of the break point considering we believe this to be a tandem duplication
left.Dup16 <- substr(genome[6], assumed.Dup16.end.position-100, assumed.Dup16.end.position-1)
right.Dup16 <- substr(genome[6], assumed.Dup16.start.position+1, assumed.Dup16.start.position+100)
# The assumed end point (15244755) works quite well. The clipped sequences before that position align to the same area
# as the FA reads (there is a short unknown sequence (AAGTAAAGTAAA) and then a sequence that aligns to pos 15222810). 
# However, I can't find any clipped reads that align at or near that point (the assumed end point was the closest I could
# get, but the clipped reads don't align at the assumed end point). 

##################
# Genotype calls #
##################

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
make.genotype.calls <- function(sample.name, method = 'new', window.buffer = 0, minimum.windows = 3, verbose = F, ploidy = 2){
	if (! method %in% c('orig', 'new'))
		stop('"method" should be one of "orig" or "new".')
	if (method == 'new')
		counts.table <- counts.list[[sample.name]]
	else
		counts.table <- counts.list_orig[[sample.name]]
	sample.name <- sub('-', '_', sample.name)
	if (verbose) cat('Sample: ', sample.name, '\n')
	# This is the vector that we are going to return after modification
	return.vector <- rep(ploidy,ncol(read.based.cyp9k1.duplications))
	names(return.vector) <- colnames(read.based.cyp9k1.duplications)
	# Manually determine genomic regions for each of the duplication types:
	Dup1.region <- Dup1.bp
	Dup2.region <- Dup2.bp
	Dup3.region <- c(15240464, Dup3.end.bp)
	Dup4.region <- Dup4.bp
	Dup5.region <- Dup5.bp
	Dup6.region <- Dup6.bp
	Dup7.region <- c(15238000, 15246300)
	# There is a deletion in Dup8, so we restrict this to the range of its actual coverage increase. 
	Dup8.region <- c(Dup8.bp[1], 15243600)
	Dup9.region <- Dup9.bp
	# Dup10 has an associated extra dup outside the CYP9K1 gene, so we restrict the range to the region without
	# that duplication.
	Dup10.region <- c(15238500, Dup10.bp[2])
	Dup11.region <- Dup11.bp
	Dup12.region <- Dup12.bp
	Dup13.region <- Dup13.bp
	Dup14.region <- c(15233807, Dup14.end.bp)
	Dup15.region <- c(15233807, Dup15.end.bp)
	Dup16.region <- c(15222810, Dup16.end.bp)
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
	pos.list[['Dup12']] <- (counts.table$Position >= Dup12.region[1]) & (counts.table$Position <= Dup12.region[2] - coverage.window)
	pos.list[['Dup13']] <- (counts.table$Position >= Dup13.region[1]) & (counts.table$Position <= Dup13.region[2] - coverage.window)
	pos.list[['Dup14']] <- (counts.table$Position >= Dup14.region[1]) & (counts.table$Position <= Dup14.region[2] - coverage.window)
	pos.list[['Dup15']] <- (counts.table$Position >= Dup15.region[1]) & (counts.table$Position <= Dup15.region[2] - coverage.window)
	pos.list[['Dup16']] <- (counts.table$Position >= Dup16.region[1]) & (counts.table$Position <= Dup16.region[2] - coverage.window)
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
					this.dup.temp.coverage <- this.dup.temp.coverage - pos.list[[this.other.dup]]*(output.rv[this.other.dup] - ploidy)
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
	sample.dups <- colnames(read.based.cyp9k1.duplications)[read.based.cyp9k1.duplications[sample.name, ] >= 1]
	sample.dups <- sample.dups[sample.dups != 'Dup0']
	remaining.dups <- sample.dups
	analysed.dups <- character()
	# If we have a Dup0, then there is something going on that we don't have a handle on so we bail
	Dup0 <- read.based.cyp9k1.duplications[sample.name, 'Dup0']
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

list.of.males <- rownames(meta)[meta$sex == 'M']
list.of.females <- rownames(meta)[meta$sex == 'F']

# Get the genotypes for all the samples
genotype.call.1 <- make.genotype.calls(sample.names[1], verbose = T)
genotype.calls <- matrix(0, nrow = length(sample.names), ncol = length(genotype.call.1), dimnames = list(sample.names, names(genotype.call.1)))
genotype.calls[sample.names[1], ] <- genotype.call.1
for (s in sample.names[2:length(sample.names)]){
	if (s %in% list.of.males)
		genotype.calls[s, ] <- make.genotype.calls(s, verbose = T, ploidy = 1)
	else
		genotype.calls[s, ] <- make.genotype.calls(s, verbose = T)
}

# Dup1
Dup1.pop.table <- table(genotype.calls[,'Dup1'], meta$population)

# Dup2
Dup2.pop.table <- table(genotype.calls[,'Dup2'], meta$population)

# Dup3:
Dup3.pop.table <- table(genotype.calls[,'Dup3'], meta$population)

# Dup4:
Dup4.pop.table <- table(genotype.calls[,'Dup4'], meta$population)
Dup4.BFS.reg.table <- table(genotype.calls[meta$population == 'BFgam','Dup4'], as.factor(as.character(meta$region[meta$population == 'BFgam'])))
Dup4.pop.table.fem <- table(genotype.calls[list.of.females, 'Dup4'], meta[list.of.females, 'population'])
HWExact(as.numeric(Dup4.pop.table.fem[,'BFgam'])) 

# Dup5:
Dup5.pop.table <- table(genotype.calls[,'Dup5'], meta$population)

# Dup6:
Dup6.pop.table <- table(genotype.calls[,'Dup6'], meta$population)
Dup6.AOM.reg.table <- table(genotype.calls[meta$population == 'AOM','Dup6'], as.factor(as.character(meta$region[meta$population == 'AOM'])))
Dup6.pop.table.fem <- table(genotype.calls[list.of.females, 'Dup6'], meta[list.of.females, 'population'])
HWExact(as.numeric(Dup6.pop.table.fem[,'AOcol'])) 

# Dup7
Dup7.pop.table <- table(genotype.calls[,'Dup7'], meta$population)

# Dup8:
Dup8.pop.table <- table(genotype.calls[,'Dup8'], meta$population)
Dup8.UG.reg.table <- table(genotype.calls[meta$population == 'UGgam','Dup8'], as.factor(as.character(meta$region[meta$population == 'UGgam'])))
Dup8.pop.table.fem <- table(genotype.calls[list.of.females, 'Dup8'], meta[list.of.females, 'population'])
HWExact(as.numeric(Dup8.pop.table.fem[,'UGgam'])) 

# Dup9:
Dup9.pop.table <- table(genotype.calls[,'Dup9'], meta$population)
Dup9.pop.table.fem <- table(genotype.calls[sub('-', '_', list.of.females), 'Dup9'], meta[list.of.females, ]$population)
HWExact(c(as.numeric(Dup9.pop.table.fem[,'CMgam']), 0)) 

# Dup10:
Dup10.pop.table <- table(genotype.calls[,'Dup10'], meta$population)

# Dup11:
Dup11.pop.table <- table(genotype.calls[,'Dup11'], meta$population)
Dup11.pop.table.fem <- table(genotype.calls[list.of.females, 'Dup11'], meta[list.of.females, 'population'])
Dup11.pop.table.mal <- table(genotype.calls[list.of.males, 'Dup11'], meta[list.of.males, 'population'])

# Dup12:
Dup12.pop.table <- table(genotype.calls[,'Dup12'], meta$population)
Dup12.GW.reg.table <- table(genotype.calls[meta$population == 'GW','Dup12'], as.factor(as.character(meta$region[meta$population == 'GW'])))
Dup12.GM.reg.table <- table(genotype.calls[meta$population == 'GM','Dup12'], as.factor(as.character(meta$region[meta$population == 'GM'])))
Dup12.pop.table.fem <- table(genotype.calls[sub('-', '_', list.of.females), 'Dup12'], meta[list.of.females, ]$population)
HWExact(as.numeric(Dup12.pop.table.fem[,'GW'])) 
HWExact(as.numeric(Dup12.pop.table.fem[,'GM'])) 

# Dup13:
Dup13.pop.table <- table(genotype.calls[,'Dup13'], meta$population)

# Dup14
Dup14.pop.table <- table(genotype.calls[,'Dup14'], meta$population)

# Dup15
Dup15.pop.table <- table(genotype.calls[,'Dup15'], meta$population)

# Dup16
Dup16.pop.table <- table(genotype.calls[,'Dup16'], meta$population)

write.table(genotype.calls, 'genotype_calls.csv', sep = '\t', col.names = NA)

save.image('CYP9K1_analysis_shrunk_data.Rdata')



