## Analyzing giant virus genomes with ViralRecall ##

First we need to download the ViralRecall repo:


> git clone https://github.com/faylward/viralrecall

and move into it:

> cd viralrecall

Now we need to download the databases used for matching:

> wget -O hmm.tar.gz https://data.lib.vt.edu/downloads/6h440s637

and unpack it:

> tar -xvzf hmm.tar.gz

For the help menu we can type:

python viralrecall.py -h


let's download the genome for Pithovirus sibericum, the largest virus discovered in terms of physical dimensions (1.5 um):
> wget -O pitho.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/916/835/GCF_000916835.1_ViralProj237323/GCF_000916835.1_ViralProj237323_genomic.fna.gz

and unzip it

> gunzip pitho.fna.gz

and then run VR on it:

> python viralrecall.py -i pitho.fna -p pitho_out -f

This will run for a while, because the HMM databases are large and the hmmsearch takes a while. Once it finishes the output files will be in the mimi_out folder. The files are:



*.faa: The proteins predicted from the input file using Prodigal

*.full_annot.tsv: A full annotation table of the predicted ORFs. This includes descriptions of the GVOG and Pfam annotations, so it can be useful if you want to look at certain annoatations in more depth.

*.vregion.annot.tsv: An annotation of only the viral regions (only present if some viral regions found)

*summary.tsv: Summary statistics for the predicted viral regions (or contig-level stats if the -c flag was used). This also includes the NCLDV marker output (marker hit: bit score)

*.pfamout: Raw output of the Pfam HMMER3 search

*.vogout: Raw output of the GVOG or VOGDB HMMER3 search

*.markerout: Raw output of the NCLDV marker gene HMMER3 search

Additionally, for each viral region viralrecall will print out .faa and .fna files for the proteins and nucleotide sequences for the regions found. Please be sure to use only .fna files as input.

Because we used the -f flag, there should also be a .pdf file there that provides the rolling average VR score across the input file contig(s). 





