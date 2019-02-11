# SEARCh
SEARCh: Sparse Enrichment Analysis for Regions in Chromatin

Usage: bash SEARCh_1.0.sh <experimental bedgraph>.bg [<control bedgraph>.bg | <FDR threshold>] ["norm" | "non"] ["union" | "AUC"]
	
	Output:

	<experimental bedgraph>.auc.threshold.merge.bed (Bed file of enriched regions)
	Data structure: <chr>	<start>	<stop>	<AUC>	<max signal>	<max signal region>

	Example:

	bash SEARCh_1.0.sh MPM091.bedgraph MPM099.bedgraph norm AUC
	Calls enriched regions in MPM091 using normalized IgG control track from MPM099
	
	bash SEARCh_1.0.sh MPM091.bedgraph MPM099.bedgraph non union
	Calls enriched regions in MPM091 using non-normalized IgG control track from MPM099 

	bash SEARCh_1.0.sh MPM091.bedgraph 0.01 non AUC
	Calls enriched regions in MPM091 by selecting the top 1% of regions by area under the curve (AUC)
