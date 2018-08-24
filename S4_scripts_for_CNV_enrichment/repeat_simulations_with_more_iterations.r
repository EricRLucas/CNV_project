load('CNV_contingencies_alt.Rdata')

# Simulate K cnv sets and record the number of cnvs that have genes, and that have detox genes. 
simulation.one.output.2 <- matrix(0,0,3, dimnames = list(c(), c('all.cnvs', 'gene.cnvs', 'detox.cnvs')))
K2 = 10000
cat('\nSimulation type 1:\n')
for (k in 1:K2){
	this.simulation <- simulate.cnvs(good.CNV.ranges.allchrom, all.chrom.positions, eligible.gene.position.indices.allchrom, detox.genes)
	simulation.one.output.2 <- rbind(simulation.one.output.2, c(nrow(this.simulation), sum(!is.na(this.simulation$Features)), sum(this.simulation$has.detox)))
	cat(k, ' ', sep = '')
}
cat('\n')
sim.above.gene.th <- simulation.one.output.2[,'gene.cnvs'] >= nrow(CNV.encompass.table.allchrom)
num.above.gene.th <- sum(sim.above.gene.th)
prop.above.gene.th <- 100 * num.above.gene.th / K2
sim.above.detox.th <- (simulation.one.output.2[,'detox.cnvs'] / simulation.one.output.2[,'gene.cnvs']) > (sum(CNV.encompass.table.allchrom$has.detox) / nrow(CNV.encompass.table.allchrom))
num.above.detox.th <- sum(sim.above.detox.th)
prop.above.detox.th <- 100 * num.above.detox.th / K2
cat('Out of ', K2, ' simulations, ', num.above.gene.th, ' (', prop.above.gene.th, '%) had at least as many cnvs containing genes as we observed in the real data. ', sep = '')
cat('In ', num.above.detox.th, ' (', prop.above.detox.th, '%) of simulations, the proportion of gene-containing cnvs that contained detox genes was at least as high as observed in the real data.\n', sep = '')


# Randomise the genes contained in CNVs and record that number that are detox genes
simulation.gene.output.2 <- matrix(0,0,3, dimnames = list(c(), c('all.clusters', 'detox.clusters', 'detox.genes')))
K2.simgene = 10000
cat('\nSimulation by gene-containing CNVs:\n')
for (k in 1:K2.simgene){
	this.simulation <- simulate.by.gene(CNV.encompass.list.allchrom, all.eligible.genes)
	simulation.gene.output.2 <- rbind(simulation.gene.output.2, c(this.simulation$all.num, this.simulation$detox.num, this.simulation$detox.gene.num))
	cat(k, ' ', sep = '')
}
cat('\n')
sim.bygene.above.th <- simulation.gene.output.2[,'detox.clusters'] >= sum(CNV.encompass.table.allchrom$has.detox)
num.bygene.above.th <- sum(sim.bygene.above.th)
prop.bygene.above.th <- 100 * num.bygene.above.th / K2.simgene
cat('Out of ', K2.simgene, ' simulations, ', num.bygene.above.th, ' (', prop.bygene.above.th, '%) had at least as many detox clusters as we observed in the real data.\n', sep = '')

save.image('repeat_simulations_with_more_iterations.Rdata')

