This directory contains test files for SEACR. The bedgraph files are as follows:

    CTCF_DE_chr1_100Mb.bedgraph.txt
    IgG_DE_chr1_100Mb.bedgraph.txt

Use SEACR to call peaks from the test files as follows:

    bash SEACR_1.3.sh CTCF_DE_chr1_100Mb.bedgraph.txt IgG_DE_chr1_100Mb.bedgraph.txt norm stringent CTCF_DE_chr1_100Mb
    bash SEACR_1.3.sh CTCF_DE_chr1_100Mb.bedgraph.txt IgG_DE_chr1_100Mb.bedgraph.txt norm relaxed CTCF_DE_chr1_100Mb

Resulting output files are as follows:

    CTCF_DE_chr1_100Mb.stringent.bed.txt
    CTCF_DE_chr1_100Mb.relaxed.bed.txt
    
The stringent output file should contain 762 lines, and the relaxed output file should contain 1146 lines. The first 10 lines of each file (output of head) should look like this:

    head CTCF_DE_chr1_100Mb.stringent.bed 
    chr1	804242	807405	1727.68	2.93081	chr1:805088-805096
    chr1	918284	921474	2315.98	3.57416	chr1:919502-919511
    chr1	936156	938790	2275.09	5.71866	chr1:937431-937432
    chr1	1307304	1308549	1407	3.35971	chr1:1307821-1307824
    chr1	1874935	1876392	1399.14	2.78784	chr1:1875346-1875376
    chr1	1890657	1892029	1349.53	4.78937	chr1:1891560-1891592
    chr1	1976474	1978693	2039.27	3.00229	chr1:1978185-1978233
    chr1	2313062	2314752	2199.11	2.93081	chr1:2313623-2313636
    chr1	2345058	2346597	1712.52	4.64641	chr1:2345791-2345935
    chr1	2477505	2480418	2908.01	3.28823	chr1:2479875-2479891

    head CTCF_DE_chr1_100Mb.relaxed.bed 
    chr1	236895	239041	879.243	0.857798	chr1:238791-238816
    chr1	619025	620892	1100.91	2.1445	chr1:620165-620171
    chr1	713215	714910	1225.86	1.85856	chr1:714127-714173
    chr1	804242	807405	1727.68	2.93081	chr1:805088-805096
    chr1	918284	921474	2315.98	3.57416	chr1:919502-919511
    chr1	936156	938790	2275.09	5.71866	chr1:937431-937432
    chr1	967182	968697	850.292	1.78708	chr1:968496-968497
    chr1	1051210	1052279	833.851	1.93005	chr1:1051566-1051604
    chr1	1226972	1228602	941.076	2.07301	chr1:1227417-1227425
    chr1	1307304	1308549	1407	3.35971	chr1:1307821-1307824
    
