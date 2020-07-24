
## Download the data

Let's start by comparing two genomes that are closely related: Escherichai coli strain K12 and Escherichia coli O157. You can probably guess that since these two genomes are from different strains in the same species that they will be pretty similar. 

promer k12.fna sent.fna

wget -O k12.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

wget -O o157.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz


mummer k12.fna o157.fna  > k12_vs_o157.mums
mummerplot --color --postscript k12_vs_o157.mums 


Now let's see what the pattern looks like as we choose genomes that are not that closely related. 

wget -O sent.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
gunzip *.gz

mummer k12.fna sent.fna  > k12_vs_sent.mums
mummerplot --color --postscript k12_vs_sent.mums -p k12_vs_sent

promer k12.fna sent.fna -p ecoli_vs_salmonella 
mummerplot --color --postscript out.delta -p ecoli_vs_salmonella

module load gnuplot/5.0.0


