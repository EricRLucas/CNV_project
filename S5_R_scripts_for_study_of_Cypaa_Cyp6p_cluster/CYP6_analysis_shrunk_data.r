library(stringr)
library(Biostrings)
library(HardyWeinberg)

load('CNV_analysis_phase2_for_CYP6.Rdata')

# Get the coordinates of all the genes in the cluster
cyp6p1.coord <- c(28499251, 28499251 + 1649)
cyp6p2.coord <- c(28501033, 28501033 + 1877)
cyp6p3.coord <- c(28491415, 28491415 + 1726)
cyp6p4.coord <- c(28497087, 28497087 + 1587)
cyp6p5.coord <- c(28494017, 28494017 + 1628)
cyp6aa1.coord <- c(28480576, 28480576 + 2061)
cyp6aa2.coord <- c(28483301, 28483301 + 1620)
cyp6ad1.coord <- c(28504248, 28504248 + 1568)
coeae60.coord <- c(28485262, 28485262 + 1818)
agap002864.coord <- c(28487640, 28487640 + 1452)

cyp6.region.coord <- c(28460000, 28570000) 
cyp6.region.plotting.indices <- c(min(which(counts.list[[1]]$Position >= cyp6.region.coord[1])), max(which(counts.list[[1]]$Position <= cyp6.region.coord[2])))

# Identify the samples with a detected duplication in the various genes in the region
encompass.duplication.detected.in.cyp6aa1 <- encompass.output.by.gene[['AGAP002862']]
encompass.duplication.detected.in.cyp6aa2 <- encompass.output.by.gene[['AGAP013128']]
encompass.duplication.detected.in.cyp6p1 <- encompass.output.by.gene[['AGAP002868']]
encompass.duplication.detected.in.cyp6p2 <- encompass.output.by.gene[['AGAP002869']]
encompass.duplication.detected.in.cyp6p3 <- encompass.output.by.gene[['AGAP002865']]
encompass.duplication.detected.in.cyp6p4 <- encompass.output.by.gene[['AGAP002867']]
encompass.duplication.detected.in.cyp6p5 <- encompass.output.by.gene[['AGAP002866']]
encompass.duplication.detected.in.cyp6.region <- unique(c(encompass.duplication.detected.in.cyp6aa1,
                                                           encompass.duplication.detected.in.cyp6aa2,
                                                           encompass.duplication.detected.in.cyp6p1,
                                                           encompass.duplication.detected.in.cyp6p2,
                                                           encompass.duplication.detected.in.cyp6p3,
                                                           encompass.duplication.detected.in.cyp6p4,
                                                           encompass.duplication.detected.in.cyp6p5))

# Load table containing geographical id
meta <- read.table('/home/eric/Liverpool/phase2.AR1/samples/samples.meta.txt', header = T, row.names = 1, sep = '\t', quote = "", comment.char = "")
rownames(meta) <- sub('-', '_', rownames(meta))
# Create booleans vector that tell us whether each individual in the samples table has a duplication in cyp6aa1
# or cyp6p3
hasencomp.aa1 <- rownames(meta) %in% encompass.duplication.detected.in.cyp6aa1
hasencomp.p3 <- rownames(meta) %in% encompass.duplication.detected.in.cyp6p3
# add these data to the meta table
meta$hasencomp.aa1 = hasencomp.aa1
meta$hasencomp.p3 = hasencomp.p3

# Now, within each population, we want to know the proportion of individuals carrying each of those duplications
propencomp.aa1 <- tapply(meta$hasencomp.aa1, meta$population, mean)
propencomp.p3 <- tapply(meta$hasencomp.p3, meta$population, mean)
propdup.for.plot.aa1 <- rbind(propencomp.aa1, 1 - propencomp.aa1)
propdup.for.plot.p3 <- rbind(propencomp.p3,  1 - propencomp.p3)
pop.sizes <- tapply(meta$population, meta$population, length)

x11(width = 15)
par(mfrow = c(2,1), mar = c(2,4,4,8), oma = c(3,0,0,0), xpd = F)
barplot(propdup.for.plot.aa1, legend.text = c('Duplication', 'No duplication'), ylab = 'Proportion', names.arg = rep('', ncol(propdup.for.plot.aa1)), main = 'CYP6AA1', args.legend = list(x=22.2, y=0.9))
barplot(propdup.for.plot.p3, ylab = 'Proportion', names.arg = rep('', ncol(propdup.for.plot.p3)),  main = 'CYP6P3')
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
SSFA.folder <- paste('/home/eric/Liverpool/CNV_/SSFA/v3_', chrom, '/phase2/CYP6AA1_region', sep = '')
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

# Get the faceaway reads in the cyp6 region
FA.cyp6.list <- list()
for (i in 1:length(FA.list)){
	this.list <- FA.list[[i]]
	this.name <- names(FA.list)[i]
	# Use the following line if you want reads that start OR end in the cyp6 region
	which.in.cyp6.FA <- ((this.list$Position >= cyp6.region.coord[1]) & (this.list$Position <= cyp6.region.coord[2])) | ((this.list$Mate.position >= cyp6.region.coord[1]) & (this.list$Mate.position <= cyp6.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the cyp6 region
#	which.in.cyp6.FA <- ((this.list$Position >= cyp6.region.coord[1]) & (this.list$Mate.position <= cyp6.region.coord[2]))
	FA.cyp6.list[[this.name]] <- this.list[which.in.cyp6.FA, ]
}
SS.cyp6.list <- list()
for (i in 1:length(SS.list)){
	this.list <- SS.list[[i]]
	this.name <- names(SS.list)[i]
	# Use the following line if you want reads that start OR end in the cyp6 region
	which.in.cyp6.SS <- ((this.list$Position >= cyp6.region.coord[1]) & (this.list$Position <= cyp6.region.coord[2])) | ((this.list$Mate.position >= cyp6.region.coord[1]) & (this.list$Mate.position <= cyp6.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the cyp6 region
#	which.in.cyp6.SS <- ((this.list$Position >= cyp6.region.coord[1]) & (this.list$Mate.position <= cyp6.region.coord[2]))
	SS.cyp6.list[[this.name]] <- this.list[which.in.cyp6.SS, ]
}
FM.cyp6.list <- list()
for (i in 1:length(FM.list)){
	this.list <- FM.list[[i]]
	this.name <- names(FM.list)[i]
	# Use the following line if you want reads that start OR end in the cyp6 region
	which.in.cyp6.FM <- ((this.list$Position >= cyp6.region.coord[1]) & (this.list$Position <= cyp6.region.coord[2])) | ((this.list$Mate.position >= cyp6.region.coord[1]) & (this.list$Mate.position <= cyp6.region.coord[2]))
	# Use the following line instead if you want reads that both start AND end in the cyp6 region
#	which.in.cyp6.FM <- ((this.list$Position >= cyp6.region.coord[1]) & (this.list$Mate.position <= cyp6.region.coord[2]))
	FM.cyp6.list[[this.name]] <- this.list[which.in.cyp6.FM, ]
}
crosschrom.cyp6.list <- list()
for (i in 1:length(crosschrom.list)){
	this.list <- crosschrom.list[[i]]
	this.name <- names(crosschrom.list)[i]
	which.in.cyp6.crosschrom <- (this.list$Position >= cyp6.region.coord[1]) & (this.list$Position <= cyp6.region.coord[2])
	crosschrom.cyp6.list[[this.name]] <- this.list[which.in.cyp6.crosschrom, ]
}

# Load up the results of the breakpoint detection 
cat('Loading CYP6 region breakpoints data\n')
bp.folder <- paste('/home/eric/Liverpool/CNV_/breakpoints/' , chrom, '/phase2/CYP6_region', sep = '')
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
clipping.start.point.cyp6.list <- list()
for (i in 1:length(clipping.start.point.list)){
	this.list <- clipping.start.point.list[[i]]
	this.name <- names(clipping.start.point.list)[i]
	which.in.cyp6.clipping.start.point <- (this.list$Position >= cyp6.region.coord[1]) & (this.list$Position <= cyp6.region.coord[2])
	this.initial.list <- this.list[which.in.cyp6.clipping.start.point, , drop = F]
	# The following for lines allow breakpoints to be kept only if their frequency exceeds a certain threshold. They 
	# are effectively pointless unless you set the threshold below above 0.
	reads.per.start.point <- tapply(this.initial.list$Position, as.factor(this.initial.list$Position), length)
	positions.to.keep <- this.initial.list$Position %in% names(reads.per.start.point)[reads.per.start.point > 0]
	clipping.start.point.cyp6.list[[this.name]] <- this.initial.list[positions.to.keep, ]
}
clipping.end.point.cyp6.list <- list()
for (i in 1:length(clipping.end.point.list)){
	this.list <- clipping.end.point.list[[i]]
	this.name <- names(clipping.end.point.list)[i]
	which.in.cyp6.clipping.end.point <- (this.list$Position >= cyp6.region.coord[1]) & (this.list$Position <= cyp6.region.coord[2])
	this.initial.list <- this.list[which.in.cyp6.clipping.end.point, , drop = F]
	reads.per.end.point <- tapply(this.initial.list$Position, as.factor(this.initial.list$Position), length)
	positions.to.keep <- this.initial.list$Position %in% names(reads.per.end.point)[reads.per.end.point > 0]
	clipping.end.point.cyp6.list[[this.name]] <- this.initial.list[positions.to.keep, ]
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
mapq0.proportions <- subset(mapq0.proportions, (Position >= cyp6.region.coord[1]) & (Position <= cyp6.region.coord[2]))
total.coverage <- mapq0.proportions$'Count mapq > 0' + mapq0.proportions$'Count mapq = 0'
mapq0.proportions$mapq0.prop <- mapq0.proportions$'Count mapq = 0' / total.coverage
# In some cases, there is no coverage, so we get infinite values. Let's deal with those
mapq0.proportions$mapq0.prop[total.coverage == 0] <- 0

# Create a list of colours that will be used for plotting the duplication types
duplication.colours <- c(rgb(0.6509804, 0.8078431, 0.8901961), rgb(0.1215686, 0.4705882, 0.7058824), rgb(0.6980392, 0.8745098, 0.5411765), rgb(0.2, 0.627451, 0.172549), rgb(0.9843137, 0.6039216, 0.6), rgb(0.8901961, 0.1019608, 0.1098039), rgb(0.9921569, 0.7490196, 0.4352941), rgb(1, 0.4980392, 0), rgb(0.7921569, 0.6980392, 0.8392157), rgb(0.4156863, 0.2392157, 0.6039216), rgb(1,0,1,0.7), rgb(0.6941176, 0.3490196, 0.1568627), rgb(0.5,0.5,0,0.8), rgb(0,0.5,0.5,0.8), rgb(0.5, 0, 0.5, 0.8))

plot.cyp6 <- function(this.sample, list.of.dups = NULL, FA = T, SS=T, FM=T, XC=T, BP=T, start.index = cyp6.region.plotting.indices[1], end.index = cyp6.region.plotting.indices[2]){
	this.sample.name <- regmatches(this.sample, regexpr('A.\\d\\d\\d\\d_C', this.sample))
	plot(counts.list[[this.sample]]$Position[start.index : end.index], counts.list[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		if (list.of.dups['Dup14'] >= 1){
			rect(Dup14.bp[1], -0.6, Dup14.bp[2], 0, col = duplication.colours[14], border = col)
			text(Dup14.bp[2], -0.4, 'Dup14', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup15'] >= 1){
			rect(Dup15.bp, -0.6, 28555300, 0, col = duplication.colours[15], border = col)
			text(28555300, -0.4, 'Dup15', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup13'] >= 1){
			rect(Dup13.bp[1], -0.6, Dup13.bp[2], 0, col = duplication.colours[13], border = col)
			text(Dup13.bp[2], -0.4, 'Dup13', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup12'] >= 1){
			rect(Dup12.bp[1], -0.6, Dup12.bp[2], 0, col = duplication.colours[12], border = col)
			text(Dup12.bp[2], -0.4, 'Dup12', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup11'] >= 1){
			rect(Dup11.bp[1], -0.6, Dup11.bp[2], 0, col = duplication.colours[11], border = col)
			text(Dup11.bp[2], -0.4, 'Dup11', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup10'] >= 1){
			rect(Dup10.bp[1], -0.6, Dup10.bp[2], 0, col = duplication.colours[10], border = col)
			text(Dup10.bp[2], -0.4, 'Dup10', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup9'] >= 1){
			rect(Dup9.bp[1], -0.6, Dup9.bp[2], 0, col = duplication.colours[9], border = col)
			text(Dup9.bp[2], -0.4, 'Dup9', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup8'] >= 1){
			rect(Dup8.bp[1], -0.6, Dup8.bp[2], 0, col = duplication.colours[8], border = col)
			text(Dup8.bp[2], -0.4, 'Dup8', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup7'] >= 1){
			rect(Dup7.bp[1], -0.6, Dup7.bp[2], 0, col = duplication.colours[7], border = col)
			text(Dup7.bp[2], -0.4, 'Dup7', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup6'] >= 1){
			rect(Dup6.bp[1], -0.6, Dup6.bp[2], 0, col = duplication.colours[6], border = col)
			text(Dup6.bp[2], -0.4, 'Dup6', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup5'] >= 1){
			rect(Dup5.bp[1], -0.6, Dup5.bp[2], 0, col = duplication.colours[5], border = col)
			text(Dup5.bp[2], -0.4, 'Dup5', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup4'] >= 1){
			rect(Dup4.bp[1], -0.6, Dup4.bp[2], 0, col = duplication.colours[4], border = col)
			text(Dup4.bp[2], -0.4, 'Dup4', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup3'] >= 1){
			rect(Dup3.bp[1], -0.6, Dup3.bp[2], 0, col = duplication.colours[3], border = col)
			text(Dup3.bp[2], -0.4, 'Dup3', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup2'] >= 1){
			rect(Dup2.bp[1], -0.6, Dup2.bp[2], 0, col = duplication.colours[2], border = col)
			text(Dup2.bp[2], -0.4, 'Dup2', cex = 1.2, adj = c(1,0))
		}
		if (list.of.dups['Dup1'] >= 1){
			rect(Dup1.bp[1], -0.6, Dup1.bp[2], 0, col = duplication.colours[8], border = col)
			text(Dup1.bp[2], -0.4, 'Dup1', cex = 1.2, adj = c(1,0))
		}
	}
	if (FA){
		if (nrow(FA.cyp6.list[[this.sample.name]]) > 0)
			add.faceaways(FA.cyp6.list[[this.sample.name]], this.col = 'blue', yrange = c(0,3))
	}
	if (SS){
		if (nrow(SS.cyp6.list[[this.sample.name]]) > 0)
			add.faceaways(SS.cyp6.list[[this.sample.name]], this.col = 'cyan', yrange = c(3,6))
	}
	if (FM){
		if (nrow(FM.cyp6.list[[this.sample.name]]) > 0)
			add.faceaways(FM.cyp6.list[[this.sample.name]], this.col = 'green', yrange = c(6,9))
	}
	if (XC){
		if (nrow(crosschrom.cyp6.list[[this.sample.name]])) 
			add.breakpoints(crosschrom.cyp6.list[[this.sample.name]]$Position, this.col = 'red', yrange = c(9,12))
	}
	if (BP){
		if (nrow(clipping.end.point.cyp6.list[[this.sample.name]]))
			add.breakpoints(clipping.end.point.cyp6.list[[this.sample.name]][, 'Position'], this.col = 'brown', yrange = c(0,12))
		if (nrow(clipping.start.point.cyp6.list[[this.sample.name]]))
			add.breakpoints(clipping.start.point.cyp6.list[[this.sample.name]][, 'Position'], this.col = 'pink', yrange = c(0,12))
	}
	these.cnv.states <- counts.list[[this.sample]]$CNV[start.index : end.index]
	lines(counts.list[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2)
	abline(v = cyp6aa1.coord[1])
	abline(v = cyp6aa1.coord[2])
	text(mean(cyp6aa1.coord), 9, 'CYP6AA1', srt=90, adj = 0, col = 'black')
	abline(v = cyp6aa2.coord[1], col = 'purple')
	abline(v = cyp6aa2.coord[2], col = 'purple')
	text(mean(cyp6aa2.coord), 9, 'CYP6AA2', srt=90, adj = 0, col = 'purple')
	abline(v = coeae60.coord[1], col = 'orange')
	abline(v = coeae60.coord[2], col = 'orange')
	text(mean(coeae60.coord), 9, 'COEAE60', srt=90, adj = 0, col = 'orange')
	abline(v = agap002864.coord[1], col = 'magenta')
	abline(v = agap002864.coord[2], col = 'magenta')
	text(mean(agap002864.coord), 9, 'AGAP002864', srt=90, adj = 0, col = 'magenta')
	abline(v = cyp6p1.coord[1], col = 'brown')
	abline(v = cyp6p1.coord[2], col = 'brown')
	text(mean(cyp6p1.coord), 9, 'CYP6P1', srt=90, adj = 0, col = 'brown')
	abline(v = cyp6p2.coord[1], col = 'green')
	abline(v = cyp6p2.coord[2], col = 'green')
	text(mean(cyp6p2.coord), 9, 'CYP6P2', srt=90, adj = 0, col = 'green')
	abline(v = cyp6p3.coord[1], col = 'violet')
	abline(v = cyp6p3.coord[2], col = 'violet')
	text(mean(cyp6p3.coord), 9, 'CYP6P3', srt=90, adj = 0, col = 'violet')
	abline(v = cyp6p4.coord[1], col = 'red')
	abline(v = cyp6p4.coord[2], col = 'red')
	text(mean(cyp6p4.coord), 9, 'CYP6P4', srt=90, adj = 0, col = 'red')
	abline(v = cyp6p5.coord[1], col = 'blue')
	abline(v = cyp6p5.coord[2], col = 'blue')
	text(mean(cyp6p5.coord), 9, 'CYP6P5', srt=90, adj = 0, col = 'blue')
}


plot.all.cyp6 <- function(list.of.all, list.of.dupl = encompass.duplication.detected.in.cyp6.region, FA=T, SS=T, FM=T, XC=T, BP=T, start.index = cyp6.region.plotting.indices[1], end.index = cyp6.region.plotting.indices[2]){
	x11()
	list.of.all <- sub('-', '_', list.of.all)
	list.of.dupl <- sub('-', '_', list.of.dupl)
	for (this.sample in list.of.all){
		if (this.sample %in% list.of.dupl)
			par(bg = 'beige')
		else
			par(bg = 'white')
		plot.cyp6(this.sample, NULL, FA, SS, FM, XC, BP, start.index, end.index)
		locator(1)
	}
}

plot.all.cyp6.with.read.dups <- function(list.of.samples, matrix.of.read.dups = read.based.cyp6.duplications, FA = T, SS = T, FM = T, XC=T, BP=T, start.index = cyp6.region.plotting.indices[1], end.index = cyp6.region.plotting.indices[2]){
	x11()
	x.midpoint <- mean(c(counts.list[[1]]$Position[start.index], counts.list[[1]]$Position[end.index]))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- sub('-', '_', list.of.samples[i])
		plot.cyp6(this.sample, matrix.of.read.dups[this.sample,], FA, SS, FM, XC, BP, start.index, end.index)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

# create a function that produces an overlapping line plot of coverage for a series of individuals in a given
# region
elegant.plot <- function(sample.names, region.coords = cyp6.region.coord, smoothing = 5, maxy = 8, ...){
	plot(region.coords, c(0,maxy), type = 'n', bty = 'n', xaxt = 'n', xlab = '', ylab = 'Normalised coverage')
	axis(1, lwd = 0, mgp = c(0,0.25,2))
	mtext(paste('Position on chromosome', chrom), 1, 1.75)
	region.indices <- which(counts.list[[1]]$Position > region.coords[1] & counts.list[[1]]$Position < region.coords[2])
	abline(h = -0.1, lwd = 2)
	rect(cyp6aa1.coord[1], -100, cyp6aa1.coord[2], 100, col = 'grey80', border = 'grey80')
	rect(cyp6aa1.coord[1], -0.2, cyp6aa1.coord[2], 0, col = 'black')
	text(mean(cyp6aa1.coord), 0.1, 'CYP6AA1', srt=90, adj = 0, col = 'black', cex = 0.8)
	rect(cyp6aa2.coord[1], -100, cyp6aa2.coord[2], 100, col = 'grey60', border = 'grey60')
	rect(cyp6aa2.coord[1], -0.2, cyp6aa2.coord[2], 0, col = 'black')
	text(mean(cyp6aa2.coord), 0.1, 'CYP6AA2', srt=90, adj = 0, col = 'black', cex = 0.8)
	rect(coeae60.coord[1], -100, coeae60.coord[2], 100, col = 'grey80', border = 'grey80')
	rect(coeae60.coord[1], -0.2, coeae60.coord[2], 0, col = 'black')
	text(mean(coeae60.coord), 0.1, 'COEAE60', srt=90, adj = 0, col = 'black', cex = 0.8)
	rect(agap002864.coord[1], -100, agap002864.coord[2], 100, col = 'grey60', border = 'grey60')
	rect(agap002864.coord[1], -0.2, agap002864.coord[2], 0, col = 'black')
	text(mean(agap002864.coord), 0.1, 'AGAP002864', srt=90, adj = 0, col = 'black', cex = 0.8)
	rect(cyp6p1.coord[1], -100, cyp6p1.coord[2], 100, col = 'grey60', border = 'grey60')
	rect(cyp6p1.coord[1], -0.2, cyp6p1.coord[2], 0, col = 'black')
	text(mean(cyp6p1.coord), 0.1, 'CYP6P1', srt=90, adj = 0, col = 'black', cex = 0.8)
	rect(cyp6p2.coord[1], -100, cyp6p2.coord[2], 100, col = 'grey80', border = 'grey80')
	rect(cyp6p2.coord[1], -0.2, cyp6p2.coord[2], 0, col = 'black')
	text(mean(cyp6p2.coord), 0.1, 'CYP6P2', srt=90, adj = 0, col = 'black', cex = 0.8)
	rect(cyp6p3.coord[1], -100, cyp6p3.coord[2], 100, col = 'grey80', border = 'grey80')
	rect(cyp6p3.coord[1], -0.2, cyp6p3.coord[2], 0, col = 'black')
	text(mean(cyp6p3.coord), 0.1, 'CYP6P3', srt=90, adj = 0, col = 'black', cex = 0.8)
	rect(cyp6p4.coord[1], -100, cyp6p4.coord[2], 100, col = 'grey80', border = 'grey80')
	rect(cyp6p4.coord[1], -0.2, cyp6p4.coord[2], 0, col = 'black')
	text(mean(cyp6p4.coord), 0.1, 'CYP6P4', srt=90, adj = 0, col = 'black', cex = 0.8)
	rect(cyp6p5.coord[1], -100, cyp6p5.coord[2], 100, col = 'grey60', border = 'grey60')
	rect(cyp6p5.coord[1], -0.2, cyp6p5.coord[2], 0., col = 'black')
	text(mean(cyp6p5.coord), 0.1, 'CYP6P5', srt=90, adj = 0, col = 'black', cex = 0.8)
	rect(cyp6ad1.coord[1], -100, cyp6ad1.coord[2], 100, col = 'grey80', border = 'grey80')
	rect(cyp6ad1.coord[1], -0.2, cyp6ad1.coord[2], 0, col = 'black')
	text(mean(cyp6ad1.coord), 0.1, 'CYP6AD1', srt=90, adj = 0, col = 'black', cex = 0.8)
	for (i in 1:length(sample.names)){
		s <- sub('-', '_', sample.names[i])
		# If needed, create a table of smoothed data here
		if (smoothing > 1){
			if (smoothing > length(region.indices))
				stop('Fail. Smoothin window size is larger than the number of points to smooth over.')
			these.counts <- counts.list[[s]][0,]
			for (j in 1:(length(region.indices) - smoothing + 1))
				these.counts[j,] <- apply(counts.list[[s]][region.indices[j:(j + smoothing - 1)], ], 2, mean)
		}
		else{
			these.counts <- counts.list[[s]][region.indices, ]
		}
		lines(these.counts$Position, these.counts$Normalised_coverage, ...)
	}
}

# Use the following definitions for the different duplications (numbers represent start range (1st row) and end range
# (2nd row) of the read pairs). These values were determined after visual inspection of the data. 
Dup1.FA = matrix(c(28480150, 28483300, 28480450, 28483600), 2, 2)
Dup1.bp = c(28480189, 28483475)
Dup1.bp.seq = c('CGTAG', 'AATTG')
#
Dup2.FA = matrix(c(28493450, 28497000, 28493750, 28497300), 2, 2)
Dup2.bp = c(28493547, 28497279)
Dup2.bp.seq = c('GCCGC','TTTAA')
#
Dup3.FA = matrix(c(28479350, 28483100, 28479650, 28483400), 2, 2)
Dup3.bp = c(28479407, 28483372)
Dup3.bp.seq = c('GCTTA', 'CAAAG')
#
Dup4.FA = matrix(c(28478850, 28482750, 28479150, 28483050), 2, 2)
Dup4.bp = c(28478925, 28483069)
Dup4.bp.seq = c('TACTT', 'CATGT')
#
Dup5.FA = matrix(c(28480300, 28484200, 28480600, 28484500), 2, 2)
Dup5.bp = c(28480372, 28484518)
Dup5.bp.seq = c('AAGAG', 'ACAAA')
#
Dup6.FA = matrix(c(28478150, 28483850, 28478450, 28484150), 2, 2)
Dup6.bp = c(28478272, 28484157)
Dup6.bp.seq = c('ATCAC', 'CTAGA')
#
Dup7.SS = matrix(c(28478000, 28486000, 28478300, 28486300), 2, 2)
Dup7.bp = c(28478057, 28486036)
Dup7.bp.seq = c('AGAGC','TTTTT')
#
Dup8.FA = matrix(c(28475900, 28484700, 28476200, 28485000), 2, 2)
Dup8.bp = c(28475996, 28485005)
Dup8.bp.seq = c('AGCGA', 'CAAAT')
#
Dup9.FA = matrix(c(28479100, 28491200, 28479400, 28491500), 2, 2)
Dup9.bp = c(28479181, 28491431)
Dup9.bp.seq = c('TGTTC', 'TGTGG')
#
Dup10.FA = matrix(c(28477800, 28490850, 28478100, 28491150), 2, 2)
Dup10.bp = c(28477889, 28491215)
Dup10.bp.seq = c('TGTAG','AACTT')
#
Dup11.FA = matrix(c(28487450, 28517800, 28487750, 28518100), 2, 2)
Dup11.bp = c(28487546, 28518123)
Dup11.bp.seq = c('AACAC', 'TTATC')
#
Dup12.FA = matrix(c(28474450, 28519650, 28474750, 28519950), 2, 2)
Dup12.bp = c(28474576, 28520016)
Dup12.bp.seq = c('CCGAC', 'ACGGT')
#
Dup13.FA = matrix(c(28472650, 28522350, 28472950, 28522650), 2, 2)
Dup13.bp = c(28472728, 28522671)
Dup13.bp.seq = c('ACCGC', 'AGCTG')
#
Dup14.FA = matrix(c(28473800, 28563200, 28474100, 28563500), 2, 2)
Dup14.bp = c(28473874, 28563596)
Dup14.bp.seq = c('CCCAC', 'AGTTG')
#
Dup15.FA = matrix(c(28465600, 55958800, 28465900, 55959100), 2, 2)
Dup15.bp = 28465673
Dup15.bp.seq = 'CAGCC'

# We will store the duplication type information in a matrix
read.based.cyp6.duplications <- matrix(0, length(counts.list), 16, dimnames = list(names(counts.list), c('Dup0', 'Dup1', 'Dup2', 'Dup3', 'Dup4', 'Dup5', 'Dup6', 'Dup7', 'Dup8', 'Dup9', 'Dup10', 'Dup11', 'Dup12', 'Dup13', 'Dup14', 'Dup15')))
# Now go through each individual and determine which duplications they have based on their diagnostic reads. 
for (this.sample in names(counts.list)){
	# get the FA and SS reads encompassed in the general region
	these.FA <- FA.cyp6.list[[this.sample]]
	these.SS <- SS.cyp6.list[[this.sample]]
	these.CSP <- clipping.start.point.cyp6.list[[this.sample]]
	these.CEP <- clipping.end.point.cyp6.list[[this.sample]]
	# Count the number of supporting reads
	num.Dup1.fa.reads <- sum((these.FA$Position > Dup1.FA[1,1]) & these.FA$Position < Dup1.FA[1,2]
	                       & (these.FA$Mate.position > Dup1.FA[2,1]) & (these.FA$Mate.position < Dup1.FA[2,2]))
	Dup1.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup1.bp[1]]
	num.Dup1.start.reads <- sum(substr(reverse(Dup1.CEP.seq), 1, 5) == Dup1.bp.seq[1])
	Dup1.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup1.bp[2]]
	num.Dup1.end.reads <- sum(substr(Dup1.CSP.seq, 1, 5) == Dup1.bp.seq[2])
	num.Dup1.reads <- num.Dup1.fa.reads + num.Dup1.start.reads + num.Dup1.end.reads
	#
	num.Dup2.fa.reads <- sum((these.FA$Position > Dup2.FA[1,1]) & these.FA$Position < Dup2.FA[1,2]
	                       & (these.FA$Mate.position > Dup2.FA[2,1]) & (these.FA$Mate.position < Dup2.FA[2,2]))
	Dup2.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup2.bp[1]]
	num.Dup2.start.reads <- sum(substr(reverse(Dup2.CEP.seq), 1, 5) == Dup2.bp.seq[1])
	Dup2.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup2.bp[2]]
	num.Dup2.end.reads <- sum(substr(Dup2.CSP.seq, 1, 5) == Dup2.bp.seq[2])
	num.Dup2.reads <- num.Dup2.fa.reads + num.Dup2.start.reads + num.Dup2.end.reads
	#
	num.Dup3.fa.reads <- sum((these.FA$Position > Dup3.FA[1,1]) & these.FA$Position < Dup3.FA[1,2]
	                       & (these.FA$Mate.position > Dup3.FA[2,1]) & (these.FA$Mate.position < Dup3.FA[2,2]))
	Dup3.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup3.bp[1]]
	num.Dup3.start.reads <- sum(substr(reverse(Dup3.CEP.seq), 1, 5) == Dup3.bp.seq[1])
	Dup3.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup3.bp[2]]
	num.Dup3.end.reads <- sum(substr(Dup3.CSP.seq, 1, 5) == Dup3.bp.seq[2])
	num.Dup3.reads <- num.Dup3.fa.reads + num.Dup3.start.reads + num.Dup3.end.reads
	#
	num.Dup4.fa.reads <- sum((these.FA$Position > Dup4.FA[1,1]) & these.FA$Position < Dup4.FA[1,2]
	                      & (these.FA$Mate.position > Dup4.FA[2,1]) & (these.FA$Mate.position < Dup4.FA[2,2]))
	Dup4.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup4.bp[1]]
	num.Dup4.start.reads <- sum(substr(reverse(Dup4.CEP.seq), 1, 5) == Dup4.bp.seq[1])
	Dup4.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup4.bp[2]]
	num.Dup4.end.reads <- sum(substr(Dup4.CSP.seq, 1, 5) == Dup4.bp.seq[2])
	num.Dup4.reads <- num.Dup4.fa.reads + num.Dup4.start.reads + num.Dup4.end.reads
	#
	num.Dup5.fa.reads <- sum((these.FA$Position > Dup5.FA[1,1]) & these.FA$Position < Dup5.FA[1,2]
	                       & (these.FA$Mate.position > Dup5.FA[2,1]) & (these.FA$Mate.position < Dup5.FA[2,2]))
	Dup5.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup5.bp[1]]
	num.Dup5.start.reads <- sum(substr(reverse(Dup5.CEP.seq), 1, 5) == Dup5.bp.seq[1])
	Dup5.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup5.bp[2]]
	num.Dup5.end.reads <- sum(substr(Dup5.CSP.seq, 1, 5) == Dup5.bp.seq[2])
	num.Dup5.reads <- num.Dup5.fa.reads + num.Dup5.start.reads + num.Dup5.end.reads
	#
	num.Dup6.fa.reads <- sum((these.FA$Position > Dup6.FA[1,1]) & these.FA$Position < Dup6.FA[1,2]
	                       & (these.FA$Mate.position > Dup6.FA[2,1]) & (these.FA$Mate.position < Dup6.FA[2,2]))
	Dup6.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup6.bp[1]]
	num.Dup6.start.reads <- sum(substr(reverse(Dup6.CEP.seq), 1, 5) == Dup6.bp.seq[1])
	Dup6.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup6.bp[2]]
	num.Dup6.end.reads <- sum(substr(Dup6.CSP.seq, 1, 5) == Dup6.bp.seq[2])
	num.Dup6.reads <- num.Dup6.fa.reads + num.Dup6.start.reads + num.Dup6.end.reads
	#
	num.Dup7.fa.reads <- sum((these.SS$Position > Dup7.SS[1,1]) & these.SS$Position < Dup7.SS[1,2]
	                       & (these.SS$Mate.position > Dup7.SS[2,1]) & (these.SS$Mate.position < Dup7.SS[2,2]))
	Dup7.CEP.seq.1 <- these.CEP$Clipped_sequence[these.CEP$Position == Dup7.bp[1]]
	num.Dup7.start.reads <- sum(substr(reverse(Dup7.CEP.seq.1), 1, 5) == Dup7.bp.seq[1])
	Dup7.CEP.seq.2 <- these.CEP$Clipped_sequence[these.CEP$Position == Dup7.bp[2]]
	num.Dup7.end.reads <- sum(substr(reverse(Dup7.CEP.seq.2), 1, 5) == Dup7.bp.seq[2])
	num.Dup7.reads <- num.Dup7.fa.reads + num.Dup7.start.reads + num.Dup7.end.reads
	#
	num.Dup8.fa.reads <- sum((these.FA$Position > Dup8.FA[1,1]) & these.FA$Position < Dup8.FA[1,2]
	                       & (these.FA$Mate.position > Dup8.FA[2,1]) & (these.FA$Mate.position < Dup8.FA[2,2]))
	Dup8.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup8.bp[1]]
	num.Dup8.start.reads <- sum(substr(reverse(Dup8.CEP.seq), 1, 5) == Dup8.bp.seq[1])
	Dup8.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup8.bp[2]]
	num.Dup8.end.reads <- sum(substr(Dup8.CSP.seq, 1, 5) == Dup8.bp.seq[2])
	num.Dup8.reads <- num.Dup8.fa.reads + num.Dup8.start.reads + num.Dup8.end.reads
	#
	num.Dup9.fa.reads <- sum((these.FA$Position > Dup9.FA[1,1]) & these.FA$Position < Dup9.FA[1,2]
	                       & (these.FA$Mate.position > Dup9.FA[2,1]) & (these.FA$Mate.position < Dup9.FA[2,2]))
	Dup9.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup9.bp[1]]
	num.Dup9.start.reads <- sum(substr(reverse(Dup9.CEP.seq), 1, 5) == Dup9.bp.seq[1])
	Dup9.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup9.bp[2]]
	num.Dup9.end.reads <- sum(substr(Dup9.CSP.seq, 1, 5) == Dup9.bp.seq[2])
	num.Dup9.reads <- num.Dup9.fa.reads + num.Dup9.start.reads + num.Dup9.end.reads
	#
	num.Dup10.fa.reads <- sum((these.FA$Position > Dup10.FA[1,1]) & these.FA$Position < Dup10.FA[1,2]
	                       & (these.FA$Mate.position > Dup10.FA[2,1]) & (these.FA$Mate.position < Dup10.FA[2,2]))
	Dup10.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup10.bp[1]]
	num.Dup10.start.reads <- sum(substr(reverse(Dup10.CEP.seq), 1, 5) == Dup10.bp.seq[1])
	Dup10.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup10.bp[2]]
	num.Dup10.end.reads <- sum(substr(Dup10.CSP.seq, 1, 5) == Dup10.bp.seq[2])
	num.Dup10.reads <- num.Dup10.fa.reads + num.Dup10.start.reads + num.Dup10.end.reads
	#
	num.Dup11.fa.reads <- sum((these.FA$Position > Dup11.FA[1,1]) & these.FA$Position < Dup11.FA[1,2]
	                       & (these.FA$Mate.position > Dup11.FA[2,1]) & (these.FA$Mate.position < Dup11.FA[2,2]))
	Dup11.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup11.bp[1]]
	num.Dup11.start.reads <- sum(substr(reverse(Dup11.CEP.seq), 1, 5) == Dup11.bp.seq[1])
	Dup11.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup11.bp[2]]
	num.Dup11.end.reads <- sum(substr(Dup11.CSP.seq, 1, 5) == Dup11.bp.seq[2])
	num.Dup11.reads <- num.Dup11.fa.reads + num.Dup11.start.reads + num.Dup11.end.reads
	#
	num.Dup12.fa.reads <- sum((these.FA$Position > Dup12.FA[1,1]) & these.FA$Position < Dup12.FA[1,2]
	                       & (these.FA$Mate.position > Dup12.FA[2,1]) & (these.FA$Mate.position < Dup12.FA[2,2]))
	Dup12.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup12.bp[1]]
	num.Dup12.start.reads <- sum(substr(reverse(Dup12.CEP.seq), 1, 5) == Dup12.bp.seq[1])
	Dup12.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup12.bp[2]]
	num.Dup12.end.reads <- sum(substr(Dup12.CSP.seq, 1, 5) == Dup12.bp.seq[2])
	num.Dup12.reads <- num.Dup12.fa.reads + num.Dup12.start.reads + num.Dup12.end.reads
	#
	num.Dup13.fa.reads <- sum((these.FA$Position > Dup13.FA[1,1]) & these.FA$Position < Dup13.FA[1,2]
	                       & (these.FA$Mate.position > Dup13.FA[2,1]) & (these.FA$Mate.position < Dup13.FA[2,2]))
	Dup13.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup13.bp[1]]
	num.Dup13.start.reads <- sum(substr(reverse(Dup13.CEP.seq), 1, 5) == Dup13.bp.seq[1])
	Dup13.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup13.bp[2]]
	num.Dup13.end.reads <- sum(substr(Dup13.CSP.seq, 1, 5) == Dup13.bp.seq[2])
	num.Dup13.reads <- num.Dup13.fa.reads + num.Dup13.start.reads + num.Dup13.end.reads
	#
	num.Dup14.fa.reads <- sum((these.FA$Position > Dup14.FA[1,1]) & these.FA$Position < Dup14.FA[1,2]
	                       & (these.FA$Mate.position > Dup14.FA[2,1]) & (these.FA$Mate.position < Dup14.FA[2,2]))
	Dup14.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup14.bp[1]]
	num.Dup14.start.reads <- sum(substr(reverse(Dup14.CEP.seq), 1, 5) == Dup14.bp.seq[1])
	Dup14.CSP.seq <- these.CSP$Clipped_sequence[these.CSP$Position == Dup14.bp[2]]
	num.Dup14.end.reads <- sum(substr(Dup14.CSP.seq, 1, 5) == Dup14.bp.seq[2])
	num.Dup14.reads <- num.Dup14.fa.reads + num.Dup14.start.reads + num.Dup14.end.reads
	#
	num.Dup15.fa.reads <- sum((these.FA$Position > Dup15.FA[1,1]) & these.FA$Position < Dup15.FA[1,2]
	                      & (these.FA$Mate.position > Dup15.FA[2,1]) & (these.FA$Mate.position < Dup15.FA[2,2]))
	Dup15.CEP.seq <- these.CEP$Clipped_sequence[these.CEP$Position == Dup15.bp[1]]
	num.Dup15.start.reads <- sum(substr(reverse(Dup15.CEP.seq), 1, 5) == Dup15.bp.seq[1])
	num.Dup15.reads <- num.Dup15.fa.reads + num.Dup15.start.reads 
	# For each of these, if there are at least 2 supporting reads, score it as a duplication. If there is only one
	# supporting read, score it tentatively
	read.based.cyp6.duplications[this.sample, 2:16] <- as.logical(c(num.Dup1.reads, num.Dup2.reads, num.Dup3.reads, num.Dup4.reads, num.Dup5.reads, num.Dup6.reads, num.Dup7.reads, num.Dup8.reads, num.Dup9.reads, num.Dup10.reads, num.Dup11.reads, num.Dup12.reads, num.Dup13.reads, num.Dup14.reads, num.Dup15.reads))
	if (num.Dup1.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup1'] <- 0.5
	if (num.Dup2.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup2'] <- 0.5
	if (num.Dup3.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup3'] <- 0.5
	if (num.Dup4.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup4'] <- 0.5
	if (num.Dup5.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup5'] <- 0.5
	if (num.Dup6.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup6'] <- 0.5
	if (num.Dup7.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup7'] <- 0.5
	if (num.Dup8.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup8'] <- 0.5
	if (num.Dup9.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup9'] <- 0.5
	if (num.Dup10.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup10'] <- 0.5
	if (num.Dup11.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup11'] <- 0.5
	if (num.Dup12.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup12'] <- 0.5
	if (num.Dup13.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup13'] <- 0.5
	if (num.Dup14.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup14'] <- 0.5
	if (num.Dup15.reads == 1)
		read.based.cyp6.duplications[this.sample, 'Dup15'] <- 0.5
	# If the sample has a duplication, but is none of the specific ones, we record it as Dup0
	if ((sum(read.based.cyp6.duplications[this.sample, ] >= 1) == 0) & (this.sample %in% encompass.duplication.detected.in.cyp6.region))
		read.based.cyp6.duplications[this.sample, 'Dup0'] <- 1
}

# Compared the duplications detected by coverage or reads
cov.based.duplications.names <- encompass.duplication.detected.in.cyp6.region
read.based.duplications.names <- rownames(read.based.cyp6.duplications)[as.logical(apply(read.based.cyp6.duplications[,2:16] >= 1, 1, sum))]

cov.based.negatives <- setdiff(read.based.duplications.names, cov.based.duplications.names)
read.based.negatives <- setdiff(cov.based.duplications.names, read.based.duplications.names)

# Output a table of each duplication type and the country they are found in
Dup0.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup0'] >= 1], 'population', drop = F]
Dup1.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup1'] >= 1], 'population', drop = F]
Dup2.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup2'] >= 1], 'population', drop = F]
Dup3.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup3'] >= 1], 'population', drop = F]
Dup4.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup4'] >= 1], 'population', drop = F]
Dup5.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup5'] >= 1], 'population', drop = F]
Dup6.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup6'] >= 1], 'population', drop = F]
Dup7.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup7'] >= 1], 'population', drop = F]
Dup8.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup8'] >= 1], 'population', drop = F]
Dup9.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup9'] >= 1], 'population', drop = F]
Dup10.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup10'] >= 1], 'population', drop = F]
Dup11.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup11'] >= 1], 'population', drop = F]
Dup12.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup12'] >= 1], 'population', drop = F]
Dup13.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup13'] >= 1], 'population', drop = F]
Dup14.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup14'] >= 1], 'population', drop = F]
Dup15.meta <- meta[rownames(read.based.cyp6.duplications)[read.based.cyp6.duplications[,'Dup15'] >= 1], 'population', drop = F]

# Let's output some lists of individuals carrying each duplication
write(rownames(Dup0.meta), file = 'CYP6_Dup0_samples.txt', ncolumns = 1)
write(rownames(Dup1.meta), file = 'CYP6_Dup1_samples.txt', ncolumns = 1)
write(rownames(Dup2.meta), file = 'CYP6_Dup2_samples.txt', ncolumns = 1)
write(rownames(Dup3.meta), file = 'CYP6_Dup3_samples.txt', ncolumns = 1)
write(rownames(Dup4.meta), file = 'CYP6_Dup4_samples.txt', ncolumns = 1)
write(rownames(Dup5.meta), file = 'CYP6_Dup5_samples.txt', ncolumns = 1)
write(rownames(Dup6.meta), file = 'CYP6_Dup6_samples.txt', ncolumns = 1)
write(rownames(Dup7.meta), file = 'CYP6_Dup7_samples.txt', ncolumns = 1)
write(rownames(Dup8.meta), file = 'CYP6_Dup8_samples.txt', ncolumns = 1)
write(rownames(Dup9.meta), file = 'CYP6_Dup9_samples.txt', ncolumns = 1)
write(rownames(Dup10.meta), file = 'CYP6_Dup10_samples.txt', ncolumns = 1)
write(rownames(Dup11.meta), file = 'CYP6_Dup11_samples.txt', ncolumns = 1)
write(rownames(Dup12.meta), file = 'CYP6_Dup12_samples.txt', ncolumns = 1)
write(rownames(Dup13.meta), file = 'CYP6_Dup13_samples.txt', ncolumns = 1)
write(rownames(Dup14.meta), file = 'CYP6_Dup14_samples.txt', ncolumns = 1)
write(rownames(Dup15.meta), file = 'CYP6_Dup15_samples.txt', ncolumns = 1)

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


######################################
# PRECISELY IDENTIFYING BREAK POINTS #
######################################

# Load the genome
genome <- readDNAStringSet('/home/eric/Liverpool/AR3/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa')

# Dup1:
assumed.Dup1.start.position <- 28480189
assumed.Dup1.end.position <- 28483475
all.Dup1.start.sequences <- character()
for (sample.name in Dup1.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup1.start.position, 'Clipped_sequence']
	all.Dup1.start.sequences <- c(all.Dup1.start.sequences, these.start.sequences)
}
longest.Dup1.start.sequence <- all.Dup1.start.sequences[which.max(nchar(all.Dup1.start.sequences))]
all.Dup1.end.sequences <- character()
for (sample.name in Dup1.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup1.end.position, 'Clipped_sequence']
	all.Dup1.end.sequences <- c(all.Dup1.end.sequences, these.end.sequences)
}
longest.Dup1.end.sequence <- all.Dup1.end.sequences[which.max(nchar(all.Dup1.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup1 <- substr(genome[2], assumed.Dup1.end.position-100, assumed.Dup1.end.position-1)
right.Dup1 <- substr(genome[2], assumed.Dup1.start.position+1, assumed.Dup1.start.position+100)
# The breakpoint looks like this (the AATT is present at both ends and therefore could be either side of the breakpoint)
# GAACCGCATGCCGATGC AATT AATTGTAATTTATTGCC
#  end of the dup ^      ^ start of the dup
#                 ^      ^ 
# position 28483470      position 28480194

# Dup2:
assumed.Dup2.start.position <- 28493547
assumed.Dup2.end.position <- 28497279
all.Dup2.start.sequences <- character()
for (sample.name in Dup2.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup2.start.position, 'Clipped_sequence']
	all.Dup2.start.sequences <- c(all.Dup2.start.sequences, these.start.sequences)
}
longest.Dup2.start.sequence <- all.Dup2.start.sequences[which.max(nchar(all.Dup2.start.sequences))]
all.Dup2.end.sequences <- character()
for (sample.name in Dup2.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup2.end.position, 'Clipped_sequence']
	all.Dup2.end.sequences <- c(all.Dup2.end.sequences, these.end.sequences)
}
longest.Dup2.end.sequence <- all.Dup2.end.sequences[which.max(nchar(all.Dup2.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup2 <- substr(genome[2], assumed.Dup2.end.position-100, assumed.Dup2.end.position-1)
right.Dup2 <- substr(genome[2], assumed.Dup2.start.position+1, assumed.Dup2.start.position+100)
# The breakpoint looks like this (the AA is present at both ends and therefore could be either side of the breakpoint)
# ATGCGTGGCCCCTCGCCG AA TTTAATTATCCGCGACG
#   end of the dup ^    ^ start of the dup
#                  ^    ^ 
#  position 28497276    position 28493550

# Dup3:
assumed.Dup3.start.position <- 28479407
assumed.Dup3.end.position <- 28483372
all.Dup3.start.sequences <- character()
for (sample.name in Dup3.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup3.start.position, 'Clipped_sequence']
	all.Dup3.start.sequences <- c(all.Dup3.start.sequences, these.start.sequences)
}
longest.Dup3.start.sequence <- all.Dup3.start.sequences[which.max(nchar(all.Dup3.start.sequences))]
all.Dup3.end.sequences <- character()
for (sample.name in Dup3.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup3.end.position, 'Clipped_sequence']
	all.Dup3.end.sequences <- c(all.Dup3.end.sequences, these.end.sequences)
}
longest.Dup3.end.sequence <- all.Dup3.end.sequences[which.max(nchar(all.Dup3.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup3 <- substr(genome[2], assumed.Dup3.end.position-100, assumed.Dup3.end.position-1)
right.Dup3 <- substr(genome[2], assumed.Dup3.start.position+1, assumed.Dup3.start.position+100)
# The breakpoint looks like this (the CA is present at both ends and therefore could be either side of the breakpoint)
# GATCGCTTGGACATTCG CA CAAAGCGCGGGTGAAT
#  end of the dup ^    ^ start of the dup
#                 ^    ^ 
# position 28483369    position 28479410

# Dup4:
assumed.Dup4.start.position <- 28478925
assumed.Dup4.end.position <- 28483069
all.Dup4.start.sequences <- character()
for (sample.name in Dup4.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup4.start.position, 'Clipped_sequence']
	all.Dup4.start.sequences <- c(all.Dup4.start.sequences, these.start.sequences)
}
longest.Dup4.start.sequence <- all.Dup4.start.sequences[which.max(nchar(all.Dup4.start.sequences))]
all.Dup4.end.sequences <- character()
for (sample.name in Dup4.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup4.end.position, 'Clipped_sequence']
	all.Dup4.end.sequences <- c(all.Dup4.end.sequences, these.end.sequences)
}
longest.Dup4.end.sequence <- all.Dup4.end.sequences[which.max(nchar(all.Dup4.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup4 <- substr(genome[2], assumed.Dup4.end.position-100, assumed.Dup4.end.position-1)
right.Dup4 <- substr(genome[2], assumed.Dup4.start.position+1, assumed.Dup4.start.position+100)
# The breakpoint looks like this 
# GGCCACACCCAATTTCATGTT CATGTAATTCATGTTCATGTTTCAT GTAATTGATCAATTGTGCGTAATTA
#      end of the dup ^                           ^ start of the dup
#                     ^       inserted seq        ^ 
#     position 28483068                           position 28478922

# Dup5:
assumed.Dup5.start.position <- 28480372
assumed.Dup5.end.position <- 28484518
all.Dup5.start.sequences <- character()
for (sample.name in Dup5.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup5.start.position, 'Clipped_sequence']
	all.Dup5.start.sequences <- c(all.Dup5.start.sequences, these.start.sequences)
}
longest.Dup5.start.sequence <- all.Dup5.start.sequences[which.max(nchar(all.Dup5.start.sequences))]
all.Dup5.end.sequences <- character()
for (sample.name in Dup5.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup5.end.position, 'Clipped_sequence']
	all.Dup5.end.sequences <- c(all.Dup5.end.sequences, these.end.sequences)
}
longest.Dup5.end.sequence <- all.Dup5.end.sequences[which.max(nchar(all.Dup5.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup5 <- substr(genome[2], assumed.Dup5.end.position-100, assumed.Dup5.end.position-1)
right.Dup5 <- substr(genome[2], assumed.Dup5.start.position+1, assumed.Dup5.start.position+100)
# The breakpoint looks like this 
# GTTCGCCGCCGATCGAGAA  ACAAACGCTAAACGCTGAC
#    end of the dup ^  ^ start of the dup
#                   ^  ^ 
#   position 28484517  position 28480373

# Dup6:
assumed.Dup6.start.position <- 28478272
assumed.Dup6.end.position <- 28484157
all.Dup6.start.sequences <- character()
for (sample.name in Dup6.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup6.start.position, 'Clipped_sequence']
	all.Dup6.start.sequences <- c(all.Dup6.start.sequences, these.start.sequences)
}
longest.Dup6.start.sequence <- all.Dup6.start.sequences[which.max(nchar(all.Dup6.start.sequences))]
all.Dup6.end.sequences <- character()
for (sample.name in Dup6.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup6.end.position, 'Clipped_sequence']
	all.Dup6.end.sequences <- c(all.Dup6.end.sequences, these.end.sequences)
}
longest.Dup6.end.sequence <- all.Dup6.end.sequences[which.max(nchar(all.Dup6.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup6 <- substr(genome[2], assumed.Dup6.end.position-100, assumed.Dup6.end.position-1)
right.Dup6 <- substr(genome[2], assumed.Dup6.start.position+1, assumed.Dup6.start.position+100)
# The breakpoint looks like this
# GCCGACACATCACCGGGCA      CTA      GAAACCGTTACTACGCGATG
#    end of the dup ^ inserted seq  ^ start of the dup
#                   ^               ^ 
#   position 28484156               position 28478273

# Dup7:
assumed.Dup7.start.position <- 28478057
assumed.Dup7.end.position <- 28486036
all.Dup7.start.sequences <- character()
for (sample.name in Dup7.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup7.start.position, 'Clipped_sequence']
	all.Dup7.start.sequences <- c(all.Dup7.start.sequences, these.start.sequences)
}
longest.Dup7.start.sequence <- all.Dup7.start.sequences[which.max(nchar(all.Dup7.start.sequences))]
# Since this is an inverted duplication, we look for clipping end points again
all.Dup7.end.sequences <- character()
for (sample.name in Dup7.samples){
	these.end.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup7.end.position, 'Clipped_sequence']
	all.Dup7.end.sequences <- c(all.Dup7.end.sequences, these.end.sequences)
}
longest.Dup7.end.sequence <- all.Dup7.end.sequences[which.max(nchar(all.Dup7.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint (the align on the reverse
# strand, as expected).
# The breakpoint reads look like this:
# GTTTGATTTTGCTTATCCTTTTT ATCGAC TCTCGATCAGTTCGATCGG
#        end of the dup ^        ^ seq just after end of the dup
#                       ^        ^ 
#       position 28478064        position 28486043
# (or its reverse complement). The reverse complement of the left side of that sequence (including the middle bit) 
# aligns to 28478058 and the right side (also including the middle bit) aligns to 28486037. The 'G' in the middle 
# sequence is a mismatch for the 28486043 alignment, so it's more likely to sit on the other side.
# The breakpoint looks like this (the G is present at both ends and therefore could be either side of the breakpoint)

# Dup8:
assumed.Dup8.start.position <- 28475996
assumed.Dup8.end.position <- 28485005
all.Dup8.start.sequences <- character()
for (sample.name in Dup8.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup8.start.position, 'Clipped_sequence']
	all.Dup8.start.sequences <- c(all.Dup8.start.sequences, these.start.sequences)
}
longest.Dup8.start.sequence <- all.Dup8.start.sequences[which.max(nchar(all.Dup8.start.sequences))]
all.Dup8.end.sequences <- character()
for (sample.name in Dup8.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup8.end.position, 'Clipped_sequence']
	all.Dup8.end.sequences <- c(all.Dup8.end.sequences, these.end.sequences)
}
longest.Dup8.end.sequence <- all.Dup8.end.sequences[which.max(nchar(all.Dup8.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup8 <- substr(genome[2], assumed.Dup8.end.position-100, assumed.Dup8.end.position-1)
right.Dup8 <- substr(genome[2], assumed.Dup8.start.position+1, assumed.Dup8.start.position+100)
# The breakpoint looks like this (the G is present at both ends and therefore could be either side of the breakpoint)
# AGTGCCCACCGTGAGCGA G CAAATTAACATTCAAC
#   end of the dup ^   ^ start of the dup
#                  ^   ^ 
#  position 28485003   position 28475998

# Dup9:
assumed.Dup9.start.position <- 28479181
assumed.Dup9.end.position <- 28491431
all.Dup9.start.sequences <- character()
for (sample.name in Dup9.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup9.start.position, 'Clipped_sequence']
	all.Dup9.start.sequences <- c(all.Dup9.start.sequences, these.start.sequences)
}
longest.Dup9.start.sequence <- all.Dup9.start.sequences[which.max(nchar(all.Dup9.start.sequences))]
all.Dup9.end.sequences <- character()
for (sample.name in Dup9.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup9.end.position, 'Clipped_sequence']
	all.Dup9.end.sequences <- c(all.Dup9.end.sequences, these.end.sequences)
}
longest.Dup9.end.sequence <- all.Dup9.end.sequences[which.max(nchar(all.Dup9.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup9 <- substr(genome[2], assumed.Dup9.end.position-100, assumed.Dup9.end.position-1)
right.Dup9 <- substr(genome[2], assumed.Dup9.start.position+1, assumed.Dup9.start.position+100)
# The breakpoint looks like this 
# CGTATCTTCTTTGTGCCT     TGT      GGGATGTTGAAGGAGAG
#   end of the dup ^              ^ start of the dup
#                  ^ inserted seq ^ 
#  position 28491430              position 28479182

# Dup10:
assumed.Dup10.start.position <- 28477889
assumed.Dup10.end.position <- 28491215
all.Dup10.start.sequences <- character()
for (sample.name in Dup10.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup10.start.position, 'Clipped_sequence']
	all.Dup10.start.sequences <- c(all.Dup10.start.sequences, these.start.sequences)
}
longest.Dup10.start.sequence <- all.Dup10.start.sequences[which.max(nchar(all.Dup10.start.sequences))]
all.Dup10.end.sequences <- character()
for (sample.name in Dup10.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup10.end.position, 'Clipped_sequence']
	all.Dup10.end.sequences <- c(all.Dup10.end.sequences, these.end.sequences)
}
longest.Dup10.end.sequence <- all.Dup10.end.sequences[which.max(nchar(all.Dup10.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup10 <- substr(genome[2], assumed.Dup10.end.position-100, assumed.Dup10.end.position-1)
right.Dup10 <- substr(genome[2], assumed.Dup10.start.position+1, assumed.Dup10.start.position+100)
# The breakpoint looks like this (the G is present at both ends and therefore could be either side of the breakpoint)
# GCAGATTCGCGCAGATGT G AACTTTATGATACTTCCG
#   end of the dup ^   ^ start of the dup
#                  ^   ^ 
#  position 28491213   position 28477891

# Dup11:
assumed.Dup11.start.position <- 28487546
assumed.Dup11.end.position <- 28518123
all.Dup11.start.sequences <- character()
for (sample.name in Dup11.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup11.start.position, 'Clipped_sequence']
	all.Dup11.start.sequences <- c(all.Dup11.start.sequences, these.start.sequences)
}
longest.Dup11.start.sequence <- all.Dup11.start.sequences[which.max(nchar(all.Dup11.start.sequences))]
all.Dup11.end.sequences <- character()
for (sample.name in Dup11.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup11.end.position, 'Clipped_sequence']
	all.Dup11.end.sequences <- c(all.Dup11.end.sequences, these.end.sequences)
}
longest.Dup11.end.sequence <- all.Dup11.end.sequences[which.max(nchar(all.Dup11.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup11 <- substr(genome[2], assumed.Dup11.end.position-100, assumed.Dup11.end.position-1)
right.Dup11 <- substr(genome[2], assumed.Dup11.start.position+1, assumed.Dup11.start.position+100)
# The breakpoint looks like this (the TTGCACG is present at both ends and therefore could be either side of the breakpoint)
# GTTTACACTACACACAA  TTGCACG   TTATCTATTCCACTGCCAA
#  end of the dup ^            ^ start of the dup
#                 ^            ^ 
# position 28518115            position 28487554

# Dup12:
assumed.Dup12.start.position <- 28474576
assumed.Dup12.end.position <- 28520016
all.Dup12.start.sequences <- character()
for (sample.name in Dup12.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup12.start.position, 'Clipped_sequence']
	all.Dup12.start.sequences <- c(all.Dup12.start.sequences, these.start.sequences)
}
longest.Dup12.start.sequence <- all.Dup12.start.sequences[which.max(nchar(all.Dup12.start.sequences))]
all.Dup12.end.sequences <- character()
for (sample.name in Dup12.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup12.end.position, 'Clipped_sequence']
	all.Dup12.end.sequences <- c(all.Dup12.end.sequences, these.end.sequences)
}
longest.Dup12.end.sequence <- all.Dup12.end.sequences[which.max(nchar(all.Dup12.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup12 <- substr(genome[2], assumed.Dup12.end.position-100, assumed.Dup12.end.position-1)
right.Dup12 <- substr(genome[2], assumed.Dup12.start.position+1, assumed.Dup12.start.position+100)
# The breakpoint looks like this (the A is present at both ends and therefore could be either side of the breakpoint)
# AATCATACGGGACCAGCC A ACGGTAAGCCAGCAAAA
#   end of the dup ^   ^ start of the dup
#                  ^   ^ 
#  position 28520014   position 28474578

# Dup13:
assumed.Dup13.start.position <- 28472728
assumed.Dup13.end.position <- 28522671
all.Dup13.start.sequences <- character()
for (sample.name in Dup13.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup13.start.position, 'Clipped_sequence']
	all.Dup13.start.sequences <- c(all.Dup13.start.sequences, these.start.sequences)
}
longest.Dup13.start.sequence <- all.Dup13.start.sequences[which.max(nchar(all.Dup13.start.sequences))]
all.Dup13.end.sequences <- character()
for (sample.name in Dup13.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup13.end.position, 'Clipped_sequence']
	all.Dup13.end.sequences <- c(all.Dup13.end.sequences, these.end.sequences)
}
longest.Dup13.end.sequence <- all.Dup13.end.sequences[which.max(nchar(all.Dup13.end.sequences))]
# The breakpoint looks like this (the TCCCG is present at both ends and therefore could be either side of the breakpoint)
#   CGTACAGCATCGCCA TCCCG AAGAGCTGACGGAAGAAG
#  end of the dup ^       ^ start of the dup
#                 ^       ^ 
# position 28522662       position 28472734

# Dup14:
assumed.Dup14.start.position <- 28473874
assumed.Dup14.end.position <- 28563596
all.Dup14.start.sequences <- character()
for (sample.name in Dup14.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup14.start.position, 'Clipped_sequence']
	all.Dup14.start.sequences <- c(all.Dup14.start.sequences, these.start.sequences)
}
longest.Dup14.start.sequence <- all.Dup14.start.sequences[which.max(nchar(all.Dup14.start.sequences))]
all.Dup14.end.sequences <- character()
for (sample.name in Dup14.samples){
	these.end.sequences <- clipping.start.point.cyp6.list[[sample.name]][clipping.start.point.cyp6.list[[sample.name]]$Position == assumed.Dup14.end.position, 'Clipped_sequence']
	all.Dup14.end.sequences <- c(all.Dup14.end.sequences, these.end.sequences)
}
longest.Dup14.end.sequence <- all.Dup14.end.sequences[which.max(nchar(all.Dup14.end.sequences))]
# The two sequences align to the right place at the other side of the breakpoint
left.Dup14 <- substr(genome[2], assumed.Dup14.end.position-100, assumed.Dup14.end.position-1)
right.Dup14 <- substr(genome[2], assumed.Dup14.start.position+1, assumed.Dup14.start.position+100)
# The breakpoint looks like this (the GC is present at both ends and therefore could be either side of the breakpoint)
# CAACAACTACTATTCCACCC GC AGTTGGGCTGAGGAAGCTGC
#     end of the dup ^    ^ start of the dup
#                    ^    ^ 
#    position 28563593    position 28473877

# Dup15:
# combine the clipping start points in the roughly correct region from all the samples with Dup15
Dup15.start.region <- c(28464000, 28467000)
all.in.Dup15.start.region <- clipping.end.point.cyp6.list[[1]][1,]
for (sample.name in Dup15.samples){
	clipping.end.point.pos <- clipping.end.point.cyp6.list[[sample.name]]$Position
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][(clipping.end.point.pos > Dup15.start.region[1]) & (clipping.end.point.pos < Dup15.start.region[2]), ]
	all.in.Dup15.start.region <- rbind(all.in.Dup15.start.region, these.start.sequences)
}
Dup15.end.region <- c(55957000, 55962000)
all.in.Dup15.end.region <- clipping.start.point.5595.list[[1]][1,]
for (sample.name in Dup15.samples){
	clipping.start.point.pos <- clipping.start.point.5595.list[[sample.name]]$Position
	these.end.sequences <- clipping.start.point.5595.list[[sample.name]][(clipping.start.point.pos > Dup15.end.region[1]) & (clipping.start.point.pos < Dup15.end.region[2]), ]
	all.in.Dup15.end.region <- rbind(all.in.Dup15.end.region, these.end.sequences)
}
# Lets look at the most repeated regions
tail(sort(table(all.in.Dup15.start.region$Position)))
# The following start position looks promising as there are lots of soft-clipped reads in that position and 
# it's very close to the FA reads alignment position. 
assumed.Dup15.start.position <- 28465673
# now the end positions
tail(sort(table(all.in.Dup15.end.region$Position)))
# Nothing looks very promising there.
# Let's look at clipped bases of the soft-clipped reads at the start position. 
all.Dup15.start.sequences <- character()
for (sample.name in Dup15.samples){
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][clipping.end.point.cyp6.list[[sample.name]]$Position == assumed.Dup15.start.position, 'Clipped_sequence']
	all.Dup15.start.sequences <- c(all.Dup15.start.sequences, these.start.sequences)
}
longest.Dup15.start.sequence <- all.Dup15.start.sequences[which.max(nchar(all.Dup15.start.sequences))]
# But blasting the remaining sequence against gambiae doesn't give much. 
# What if we create a table of breakpoint reads from non-Dup15 individuals, and see how many of them have those
# clipped reads
all.non.Dup15.start.region <- clipping.end.point.cyp6.list[[1]][1,]
for (sample.name in setdiff(sample.names, Dup15.samples)){
	clipping.end.point.pos <- clipping.end.point.cyp6.list[[sample.name]]$Position
	these.start.sequences <- clipping.end.point.cyp6.list[[sample.name]][(clipping.end.point.pos > Dup15.start.region[1]) & (clipping.end.point.pos < Dup15.start.region[2]), ]
	all.non.Dup15.start.region <- rbind(all.non.Dup15.start.region, these.start.sequences)
}
# Only 5 reads were soft clipped at that position in the non-Dup15 samples, and none of those had the "right" bases
# after the clipping point, so we can assume that those soft-clipped reads at least identify the dup, if not cover
# the breakpoint. 
# Those are probably the breakpoint reads, but the other side is probably a repetitive region that 
# could align in lots of place. That's why only a few of the Dup15 samples have the supporting FA reads. Indeed, 
# looking at the crosschrom reads of, eg, sample AY0069, there are some that align to the start region of Dup15, with
# the mate aligning to other chromosomes. This supports that possibility. Still, we can now make calls with the 
# breakpoint reads. 
# The breakpoint looks like this 
# ACTCCGAACGACTCCGAC  TCGCTACGGATTCAGAA
#   end of the dup ^  ^ start of the dup
#                  ^  ^ 
#   position unknown  position 28465674

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
	return.vector <- rep(2,ncol(read.based.cyp6.duplications))
	names(return.vector) <- colnames(read.based.cyp6.duplications)
	# Manually determine genomic regions for each of the duplication types:
	Dup1.region <- Dup1.bp
	Dup2.region <- Dup2.bp
	Dup3.region <- Dup3.bp
	Dup4.region <- Dup4.bp
	Dup5.region <- Dup5.bp
	Dup6.region <- Dup6.bp
	Dup7.region <- Dup7.bp
	Dup8.region <- Dup8.bp
	Dup9.region <- Dup9.bp
	Dup10.region <- Dup10.bp
	Dup11.region <- Dup11.bp
	Dup12.region <- Dup12.bp
	Dup13.region <- Dup13.bp
	Dup14.region <- Dup14.bp
	Dup15.region <- c(Dup15.bp, 28555000)
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
	sample.dups <- colnames(read.based.cyp6.duplications)[read.based.cyp6.duplications[sample.name, ] >= 1]
	sample.dups <- sample.dups[sample.dups != 'Dup0']
	remaining.dups <- sample.dups
	analysed.dups <- character()
	# If we have a Dup0, then there is something going on that we don't have a handle on so we bail
	Dup0 <- read.based.cyp6.duplications[sample.name, 'Dup0']
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
Dup1.UG.reg.table <- table(genotype.calls[meta$population == 'UGgam','Dup1'], as.factor(as.character(meta$region[meta$population == 'UGgam'])))
HWExact(as.numeric(Dup1.pop.table[,'UGgam'][1:3])) 

# Dup2:
Dup2.pop.table <- table(genotype.calls[,'Dup2'], meta$population)

# Dup3:
Dup3.pop.table <- table(genotype.calls[,'Dup3'], meta$population)

# Dup4:
Dup4.pop.table <- table(genotype.calls[,'Dup4'], meta$population)
HWExact(as.numeric(Dup4.pop.table[,'BFcol'])) 

# Dup5:
Dup5.pop.table <- table(genotype.calls[,'Dup5'], meta$population)

# Dup6:
Dup6.pop.table <- table(genotype.calls[,'Dup6'], meta$population)

# Dup7:
Dup7.pop.table <- table(genotype.calls[,'Dup7'], meta$population)
Dup7.CI.reg.table <- table(genotype.calls[meta$population == 'CIcol','Dup7'], as.factor(as.character(meta$region[meta$population == 'CIcol'])))
Dup7.BFM.reg.table <- table(genotype.calls[meta$population == 'BFcol','Dup7'], as.factor(as.character(meta$region[meta$population == 'BFcol'])))
HWExact(as.numeric(Dup7.pop.table[,'CIcol'])) 
HWExact(as.numeric(Dup7.pop.table[,'BFcol'])) 

# Dup8:
Dup8.pop.table <- table(genotype.calls[,'Dup8'], meta$population)

# Dup9:
Dup9.pop.table <- table(genotype.calls[,'Dup9'], meta$population)

# Dup10:
Dup10.pop.table <- table(genotype.calls[,'Dup10'], meta$population)
Dup10.BFM.reg.table <- table(genotype.calls[meta$population == 'BFcol','Dup10'], as.factor(as.character(meta$region[meta$population == 'BFcol'])))
HWExact(as.numeric(Dup10.pop.table[,'BFcol'])) 

# Dup11
Dup11.pop.table <- table(genotype.calls[,'Dup11'], meta$population)

# Dup12:
Dup12.pop.table <- table(genotype.calls[,'Dup12'], meta$population)

# Dup13
Dup13.pop.table <- table(genotype.calls[,'Dup13'], meta$population)

# Dup14:
Dup14.pop.table <- table(genotype.calls[,'Dup14'], meta$population)
Dup14.CI.reg.table <- table(genotype.calls[meta$population == 'CIcol','Dup14'], as.factor(as.character(meta$region[meta$population == 'CIcol'])))
HWExact(as.numeric(Dup14.pop.table[,'CIcol'][1:3])) 

# Dup15
Dup15.pop.table <- table(genotype.calls[,'Dup15'], meta$population)
Dup15.CI.reg.table <- table(genotype.calls[meta$population == 'CIcol','Dup15'], as.factor(as.character(meta$region[meta$population == 'CIcol'])))
HWExact(as.numeric(Dup15.pop.table[,'CIcol'])) 
write.table(genotype.calls, 'genotype_calls.csv', sep = '\t', col.names = NA)

save.image('CYP6_analysis_shrunk_data.Rdata')



