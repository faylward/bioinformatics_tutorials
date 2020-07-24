## Gene Prediction ##

Let's start by downloading a Staphylococcus aureus genome from NCBI:

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/585/GCF_000009585.1_ASM958v1/GCF_000009585.1_ASM958v1_genomic.fna.gz

and let's make sure to unzip it so that we can access the .fna file directly (gunzip command). 
Make sure to use "head" and "tail" as we did in the W1 tutorial to ensure that the file is in FASTA format.

gunzip  GCF_000009585.1_ASM958v1_genomic.fna.gz

and

ls -lh

