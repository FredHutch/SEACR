This directory contains test files for SEACR. The bedgraph files are as follows:

    CTCF_DE_chr1_100Mb.bedgraph.txt
    IgG_DE_chr1_100Mb.bedgraph.txt

Use SEACR to call peaks from the test files as follows:

    bash SEACR_1.4.sh -b CTCF_DE_chr1_100Mb.bedgraph.txt -c IgG_DE_chr1_100Mb.bedgraph.txt -n norm -m stringent -o CTCF_DE_chr1_100Mb
    bash SEACR_1.4.sh -b CTCF_DE_chr1_100Mb.bedgraph.txt -c IgG_DE_chr1_100Mb.bedgraph.txt -n norm -m relaxed -o CTCF_DE_chr1_100Mb

Resulting output files are as follows:

    CTCF_DE_chr1_100Mb.stringent.bed.txt
    CTCF_DE_chr1_100Mb.relaxed.bed.txt
    
The stringent output file should contain 762 lines, and the relaxed output file should contain 1146 lines. The first 10 lines of each file (output of head) should look like this:

    head CTCF_DE_chr1_100Mb.stringent.bed 
    chr1	236895	239041	879.243	0.857798	chr1:238791-238816
    chr1	368948	370425	779.953	1.7156	chr1:369527-369556
    chr1	619025	620892	1100.91	2.1445	chr1:620165-620171
    chr1	713215	714910	1225.86	1.85856	chr1:714127-714173
    chr1	804242	807405	1727.68	2.93081	chr1:805088-805096
    chr1	918284	921474	2315.98	3.57416	chr1:919502-919511
    chr1	936156	938790	2275.09	5.71866	chr1:937431-937432
    chr1	967182	968697	850.292	1.78708	chr1:968496-968497
    chr1	1051210	1052279	833.851	1.93005	chr1:1051566-1051604
    chr1	1226972	1228602	941.076	2.07301	chr1:1227417-1227425

    head CTCF_DE_chr1_100Mb.relaxed.bed 
    chr1	10013	11929	980.606	1.7156	chr1:10346-10357
    chr1	91962	93519	491.018	1.00076	chr1:92387-92396
    chr1	105013	105797	498.737	1.07225	chr1:105214-105227
    chr1	236895	239041	879.243	0.857798	chr1:238791-238816
    chr1	251258	252108	384.65	1.00076	chr1:251347-251526
    chr1	368948	370425	779.953	1.7156	chr1:369527-369556
    chr1	415445	416696	419.749	0.929282	chr1:416448-416473
    chr1	619025	620892	1100.91	2.1445	chr1:620165-620171
    chr1	713215	714910	1225.86	1.85856	chr1:714127-714173
    chr1	780230	781942	385.437	1.14373	chr1:781160-781201

    
