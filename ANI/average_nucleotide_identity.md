## Average Nucleotide Identity ##
<br>
<p> Examining the whole-genome nucleotide identity between genomes is often a useful tool for classification. 
<br>
<br>

>mamba install fastani -c bioconda

The manual is here: https://github.com/ParBLiSS/FastANI

And now we need to download some genomes

wget -O k12.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

And

wget -O shigella.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/925/GCF_000006925.2_ASM692v2/GCF_000006925.2_ASM692v2_genomic.fna.gz

And now we can run fastANI:

> fastANI -q k12.fna.gz -r shigella.fna.gz -o ani_out.txt

We can see the results in the ani_out.txt file. The columns are query genome, reference genome, ANI value, count of bidirectional fragment mappings, and total query fragments. The alignment fraction is simply the ratio of mappings and total fragments.  
In this case the ANI is 97% and the alignment fraction is  1322/1547 = 85%. So these two genomes are quite similar!

Now let's do it the long way with blastn and compare:

> gunzip *

> prodigal -i k12.fna -d k12.genes.fna

> prodigal -i shigella.fna -d shigella.genes.fna

> makeblastdb -in shigella.genes.fna -dbtype nucl

> blastn -query k12.genes.fna -db shigella.genes.fna -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 | datamash mean 3

I got 98.4%, which is pretty close to the fastANI result. 

