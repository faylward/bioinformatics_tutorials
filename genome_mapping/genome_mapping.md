## Whole genome mapping and variant calling ##

Here we will be looking at the genome of Yersinia pestis, which can be found on NCBI at this location
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/222/975/GCF_000222975.1_ASM22297v1/

Copy this URL into your browser and take a look at the files. These are publicly-available files that are made available from the National Center for Biotechnology Information (NCBI). 

Note that many of the files are in a compressed .gz format. 

.fna files are Fasta Nucleic Acid (chromosome or gene sequences)
.faa files are Fasta Amino Acid   (protein sequences)

Here we will be most interested in the gene and chromosome files.

<br/>
<br/>

### Download some data 

To get started let's download the genes, in FASTA format, for Yersinia Pestis CO92

> wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/222/975/GCF_000222975.1_ASM22297v1/GCF_000222975.1_ASM22297v1_cds_from_genomic.fna.gz

This command uses the common Unix utility "wget", which will download a file directy to the folder in which you are located. After doing this you should see the .fna.gz file in your folder. You can check this with the "ls" command. 

<br/>
<br/>