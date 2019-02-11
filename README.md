# SEARCh
SEARCh: Sparse Enrichment Analysis for Regions in Chromatin

Usage: 

	bash SEARCh_1.0.sh experimental bedgraph [control bedgraph | FDR threshold] ["norm" | "non"] ["union" | "AUC"]
	
Output:

	<experimental bedgraph>.auc.threshold.merge.bed (Bed file of enriched regions)
Data structure: 
	
	<chr>	<start>	<stop>	<AUC>	<max signal>	<max signal region>

Example:

	bash SEARCh_1.0.sh target.bedgraph IgG.bedgraph norm AUC
Calls enriched regions in target data using normalized IgG control track using AUC threshold
	
	bash SEARCh_1.0.sh target.bedgraph IgG.bedgraph non union
Calls enriched regions in target data using non-normalized IgG control track using AUC and max signal thresholds 

	bash SEARCh_1.0.sh target.bedgraph 0.01 non AUC
Calls enriched regions in MPM091 by selecting the top 1% of regions by AUC
