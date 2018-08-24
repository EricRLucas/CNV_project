# This script reduces the output from the CNV analysis to reduce its memory footprint

load('/home/eric/Liverpool/CNV_v2_bigfiles/counting_output_v4_3R/phase2/alt_fullgcnorm_nomapq_mapq002_varvar_trans000001/CNV_analysis_alt_3R.Rdata')

cyp6m2.z1.region.coord <- c(6900000, 7000000) 

	
for (sample.name in names(counts.list)){
	cat('\t', sample.name, '\n', sep = '')
	counts.list[[sample.name]] <- subset(counts.list[[sample.name]], (Position >= cyp6m2.z1.region.coord[1]) & (Position <= cyp6m2.z1.region.coord[2]))
}
	
# Now save this to an .Rdata file
save.image('CNV_analysis_phase2_for_CYP6M2-Z1.Rdata')

