# 3983
RNA-Seq for EcoEvo #2– Modularity of gene expression in butterfly evolution

This repository contains custom scripts for NBIS support project #3983 "RNA-Seq for EcoEvo #2– Modularity of gene expression in butterfly evolution".

## Contents

* `annotate_denovo_clusters_arp7.pl` and `annotate_denovo_clusters_arp7_top5.pl` annotate positive hits in results of reciprocal blast using homologue groups from [EvidentialGene: Orthology completeness of Arthropod gene assemblies](http://arthropods.eugenes.org/EvidentialGene/arthropods/Arthropod_Orthology_Completeness/). Files used are from the ARP7 release: [arp7s10f_genes.ugp.txt](http://arthropods.eugenes.org/arthropods/orthologs/ARP7/arp7bor5b/arp7s10f_genes.ugp.txt.gz) and for annotating ARP7 groups with descriptive names [arp7s10f_omclgn.consensus_def.txt](http://arthropods.eugenes.org/arthropods/orthologs/ARP7/arp7bor5b/arp7s10f_omclgn.consensus_def.txt).

 `annotate_denovo_clusters_arp7.pl` produces annotated list of ARP7 homologues using the **top hit** in the reciprocal blast

 `annotate_denovo_clusters_arp7_top5.pl` produces annotated list of ARP7 homologues using the **any of the 5 top hits** in the reciprocal blast if it belongs to the same ARP7 homologue group as the first hit

 