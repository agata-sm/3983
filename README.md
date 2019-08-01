# 3983
RNA-Seq for EcoEvo #2– Modularity of gene expression in butterfly evolution

This repository contains custom scripts for NBIS support project #3983 "RNA-Seq for EcoEvo #2– Modularity of gene expression in butterfly evolution".

## Contents (custom scripts)

* `parse_trinotate_report_ids.pl` parses report generated by `Trinotate` to obtain a list of partial gene IDs.

* `get_ids_from_fa.pl` obtains a list of full IDs (as in the original `fasta` file) based on the partial IDs from `Trinotate` report.

* `parse_reblast_arp7_2018_v2.pl` parses results of `blast` and reciprocal `blast`

* `annotate_denovo_clusters_arp7.pl` and `annotate_denovo_clusters_arp7_top5.pl` annotate positive hits in results of reciprocal `blast` using homologue groups from [EvidentialGene: Orthology completeness of Arthropod gene assemblies](http://arthropods.eugenes.org/EvidentialGene/arthropods/Arthropod_Orthology_Completeness/). Files used are from the ARP7 release: [arp7s10f_genes.ugp.txt](http://arthropods.eugenes.org/arthropods/orthologs/ARP7/arp7bor5b/arp7s10f_genes.ugp.txt.gz) and for annotating ARP7 groups with descriptive names [arp7s10f_omclgn.consensus_def.txt](http://arthropods.eugenes.org/arthropods/orthologs/ARP7/arp7bor5b/arp7s10f_omclgn.consensus_def.txt).

 `annotate_denovo_clusters_arp7.pl` produces annotated list of ARP7 homologues using the **top hit** in the reciprocal blast;

 `annotate_denovo_clusters_arp7_top5.pl` produces annotated list of ARP7 homologues using the **any of the 5 top hits** in the reciprocal blast if it belongs to the same ARP7 homologue group as the top hit.

* `integrate_homologues_counts_v0.2.1.pl` merges ARP7 homologue group (a.k.a gene) level count tables from four butterfly species.


## Workflow

The workflow for annotation of *de novo* assembled transcripts used in this study is summarised below. Details are given in the project report, section *Bioinformatics Methods*.

1. `Trinotate` pipeline, including `BLAST` search for homologues (`blastp` using ORFs predicted by `TransDecoder.Predict`).

2. Reciprocal `BLAST` (`blastp`). All ARP7 gene IDs were used. A custom perl script `get_ids_from_fa.pl` was used to obtain full IDs
from fasta file based on a list of partial IDs output by `parse_trinotate_report_ids.pl`.


3. 
