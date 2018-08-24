# This script reduces the output from the CNV analysis to reduce its memory footprint 

load('/home/eric/Liverpool/CNV_v2_bigfiles/counting_output_v4_2R/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/CNV_analysis_alt_2R.Rdata')

cyp6.region.coord <- c(28460000, 28570000) 

	
for (sample.name in names(counts.list)){
	cat('\t', sample.name, '\n', sep = '')
	counts.list[[sample.name]] <- subset(counts.list[[sample.name]], (Position >= cyp6.region.coord[1]) & (Position <= cyp6.region.coord[2]))
}
	
# Now save this to an .Rdata file
save.image('CNV_analysis_phase2_for_CYP6.Rdata')

