# Giant_viruses
Script for analyzing MGE associated with defense systems, using data obtained from DefenseFinder or PADLOC. 

This script is used to search for MGEs that span five genes upstream and downstream of each defense system (a surrounding window of ten genes) and compare them with the MobileOG-db database. This database contains proteins that mediate functions essential for MGEs: (i) integration and excision (IE) from one genetic locus to another; (ii) nucleic acid replication, recombination, or repair (RRR); (iii) transfer between organisms (T); (iv) stability, transfer, or defense elements (STD); and (v) phage-specific biological processes (P) (e.g., genome packaging or lysis and lysogeny), as well as the transcriptional regulators associated with these processes. Next, the developed script was used to quantify the aforementioned proteins in the genomes of giant viruses. Sites were separated into three categories (integrative, prophage, or conjugative) if they matched “integration,” “phage,” or “conjugation” in their main category annotation.

To obtain mobilome data, go to: https://github.com/clb21565/mobileOG-db and cite: 

Brown, C. L., Mullet, J., Hindi, F., Stoll, J. E., Gupta, S., Choi, M., Keenum, I., Vikesland, P., Pruden, A., & Zhang, L. (2022). mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements. Applied and environmental microbiology, 88(18), e0099122. https://doi.org/10.1128/aem.00991-22
