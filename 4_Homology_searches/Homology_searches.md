## Homology searches ##
<br>
<p> When analyzing genomic data we will often find ourselves in a situation where we have a genome/gene/protein of interest and we will want to know what it is similar to. Perhaps you have a specific gene you're interested in, and you want to know what organisms have a similar gene. Or perhaps you have a new genome, and you want to know what closely-related genomes have already been sequenced. A general first step in answering these questions involves doing a simple homology search to compare your sequence of interest to a database of reference sequences. That's what we'll focus on here, with an emphasis on protein-protein comparisons. 
<br>
<p> Here we will use proteins predicted from the genomes of two Prochlorococcus bacteriophage genomes. We can use the wget command, which is already available as part of the base Ubuntu command line. Wget allows us to download files from a web server directly into the folder we are working in, and we need to know the URL for the file in order to do this. The National Center for Biotechnology Information (NCBI) has many genomes and genome-related datasets that it posts for researchers to use, and I have gone through and found the appropriate URLs to use here. <br>
<p> Here are the commands:
 
<br>

>wget -O PSSM2.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/859/585/GCF_000859585.1_ViralProj15135/GCF_000859585.1_ViralProj15135_genomic.fna.gz

>wget -O PSSM3.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/907/775/GCF_000907775.1_ViralProj209210/GCF_000907775.1_ViralProj209210_genomic.fna.gz


<p> The -O flag specifies the file names that we want the downloads to be called. Without this flag wget would give the downloaded file the same names that they are given on the website, and sometimes these names can be quite long. 
 
The above commands should download one gzip file each (extension .gz). Gzip files are commonly used for compressing data on Linux systems. Before we use them here we will have to uncompress them with the gunzip command
gunzip *.gz
 
After this we should have two .faa files. To check this we can simply use the ls command:

> ls 

 <br>
### Predict proteins

 
Now that we have our genomes downloaded and unzipped, we need to predict genes and proteins (similar to what what is described in a [Previous Tutorial](https://github.com/faylward/bioinformatics_tutorials/blob/master/3_Gene_Prediction/Gene_Prediction.md))

> prodigal -i PSSM2.fna -a PSSM2.faa -d PSSM2.genes.fna -f gff -o PSSM2.gff

> prodigal -i PSSM3.fna -a PSSM3.faa -d PSSM3.genes.fna -f gff -o PSSM3.gff
 <br>
 
### Get some statistics on the number of proteins encoded in each genome

Now that we have the files downloaded and in the right format, we can get some basic stats about their format and content. FASTA files are formatted such that sequences are always preceded by a "header" line that starts with ">". This line contains information about the name of the sequence, and possibly other information. 
Let's use "seqkit stats" to get some sequence statistics for the protein files. 

> seqkit stats PSSM2.faa 

or

> seqkit stats PSSM3.faa 

### Make a reference BLAST database to search against

Now that we have some basic information about the files we can begin formatting them for BLASTP.
For any search BLAST needs one FASTA file to be specified as the query, and one to be specified as the reference. Before running BLAST the reference needs to be formatted using a command called makeblastdb, which is part of the BLAST package and should have bee installed above. Makeblastdb takes a FASTA file as input and produces several files with different extensions that can be used as reference databases. 
 
> makeblastdb -in PSSM2.faa -dbtype prot
 
the -in flag specified the FASTA file to be formatted, and the -dbtype flag specifies the molecule type (prot for protein and nucl for nucleic acid). 
 
If we check the folder after running the above command we should see a number of new files with the prefix PSSM2.faa and several new suffixes. 

### Run a BLASTP search

Now that we have one file formatted as a reference database we can run BLASTP, using the other FASTA file as the query. 
 
> blastp -query PSSM3.faa -db PSSM2.faa | head -n 100
 
The output here is quite long, so we can pipe the output into a head command so we see only the first 100 lines. You should see something like this (note not all 100 lines are shown below). 

Note that there is a lot of information in this output. The program name and information are provided, and there is information on the query and reference databases used. The alignments that were calculated are also provided- you can see the top of one in the image above, and you can scroll to inspect it when you run it yourself. Note that this information is provided for every alignment that could be calculated for each protein in the query file, so if we had put this output into a file it would be quite large. 

### Getting a simplified output

Since the output of the last step was quite extensive, we will want to find ways to simplify it. 
Here is a similar command that will provide tab-delimited output (first 10 hits shown with the head command). 

> blastp -query PSSM3.faa -db PSSM2.faa -outfmt 6 | head

A few notes on the output here:
1) The format is tab-delimited, and the columns correspond to different statistics that were calculated for individual alignments. The columns are: query protein, reference protein, % identity, alignment length, mismatches, gap opens, query start, query end, reference start, reference end, evalue, bit score
2) Every protein in the query is compared to every protein in the reference, and all alignments are reported. So a protein in the query file could conceivable have alignments to multiple proteins in the reference (indeed, we see this in the image above). Also, multiple alignments that could be found between the same proteins are also shown (for example, if only the beginning and end of the amino acid sequences align, then two distinct alignmens will be provided). 
3) Just becuase an alignment is reported does not mean it is 'real'. There are cases of 'spurious alignments' that could happend by random chance. We will need to investigate the statistics provided for each alignment to decide whether or not we think it's worth trusting. 
 
Here are some additional parameters that will ensure that very poor alignments are not reported, and that only the best alignment for each query protein are given (for simplicity). 
blastp -query PSSM3.faa -db PSSM2.faa -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -max_hsps 1 -qcov_hsp_perc 50 | head
 
Different users will prefer different e-values and other cutoffs depending on what they are trying to do afterwards, their own comfort level, etc. As a common rule-of-thumb, e-values of 1e-3 and qe-5 are pretty common. For your own analyses you will need to use your own biological insight to decide for yourself what you are willing to trust and whether the results make sense.
 
Here is a breakdown of the flags used above:
 -query: this is the input file, so the file with all of the protein sequences that we want to search

 -db: this is the database, so the file we just indexed with the makeblastdb command above. Note that makeblastdb creates multiple reference files and that only the root name needs to be given here (so if the database was called refdb, then refdb would be given here even though the index files are called refdb.pin, refdb.phr, etc.)

 -max_target_seqs: This flat specifies that we only want the best hit for each query protein. Otherwise all hits are provided.

 -outfmt: This specifies that we want the tab-delimited output format rather than the full alignment output. If you forget what the columns are you can use -outfmt 7.

 -evalue: This indicates that we want to exclude all hits with evalues above this threshold. A good value is about 0.00001, or 1e-5.

 max_hsps: HSPs are 'high-scoring segment pairs'. A query protein can make several separate alignments to a single reference, so this tells the program we want only the best-scoring alignment.

 -qcov_hsp_perc: This is the 'query coverage high-scoring sequence pair percent', or the percent of the query protein that has to form an alignment against the reference to be retained. Higher values prevent spurious alignments of only a short portion of the query to a reference.

### Filtering the BLASTP output

Since we are comparing all of the proteins encoded in two viral genomes, it would be nice to get two basic statistics:
1) How many proteins in genome A are present in genome B and vice versa.
2) Of he proteins that are present in both genomes, how similar are they overall?
 
To answer the first question, we can simply use the command in the last step and count how many hits we find overall. Using the commands described in the previous step we need to be careful to make sure we are only counting the best hit for each query protein, and only one alignment per protein pair. 
 
> blastp -query PSSM3.faa -db PSSM2.faa -outfmt 6 -max_target_seqs 1 -evalue 0.00001 -max_hsps 1 -qcov_hsp_perc 50 | wc
 
Note that instead of piping the output to the 'head' command, as we did above, now we can pipe it into a 'wc' command to count how many output lines there are. You should see something like this:

So according to this analysis, there are 119 proteins in phage PSSM3 that have hits to phage PSSM2, with the parameters we used. What percent of all proteins in PSSM3 have hits to proteins in PSSM2? And vice versa?
 
Now as an exercise try changing the parameters a bit and see how they change the output. What do you think lowering the e-value threshold will do? What is the result of changing the query coverage percent?
 
Importantly, note that the results may change if we switch the query and the reference files (why would this be?), so we will want to do the reciprocal analysis too. 

### Plotting the Results in the R console

Above I mentioned two questions we would like to answer:
1) How many proteins in genome A are present in genome B and vice versa.
2) Of he proteins that are present in both genomes, how similar are they overall?
 
We answered the first qeustion above, and for the second we can save the output of a BLASTP search to a file and then plot the results in R. 

> blastp -query PSSM3.faa -db PSSM2.faa -outfmt 6 -max_target_seqs 1 -evalue 0.00001 -max_hsps 1 -qcov_hsp_perc 50 > PSSM3_vs_PSSM2.out

Now for the next few lines of code we're going to enter the R environment. You will need to open up an R console for this, or if you are working in an RStudio environment click on the "Consol" tab on the bottom-left portion of the screen. 

First we need to navigate to the appropriate folder again (we may have done this before in the terminal but we need to do it again in R). The R command for this is "setwd", short for "set working directory".
setwd("/home/faylward/tutorials/")

Then we import the data into an R object. We can call the object whatever we want. For purposes here we'll just call it "blastout". 
The command for reading in tabular data is "read.table", and you can see we need extra flags in the command to specify the file name and that the columns are tab-delimited. 
>blastout <- read.table(file="PSSM3_vs_PSSM2.out", sep="\t")

We can look inside the blastout object we created using a command similar to the "head" command in Unix. 
>head(blastout)

We can make a simple histogram of the Percent Identity column using the "hist" command. Notice how we specify the correct column- we need to use the $ symbol followed by the correct column header. 
>hist(blastout$V3)

We can play around with various other plotting parameters to make the plot look nicer. R is really good for this. Here are a few samples:
>hist(blastout$V3, breaks=50, col="blue")

>hist(blastout$V3, breaks=50, col="blue", main="PSSM3 vs PSSM2 BLASTP % ID", xlab="% ID", ylab="number of hits")

>hist(blastout$V3, breaks=50, col="blue", main="PSSM3 vs PSSM2 BLASTP % ID", xlab="% ID", ylab="number of hits", xlim=c(0, 100))
 
![Image of Example RStudio Console](https://github.com/faylward/bioinformatics_tutorials/blob/master/4_Homology_searches/image10.png)
 

### Calculating the Average Amino Acid Identity
For closely related genomes many scientists prefer using average nucleic acid identity (ANI) instead, but for distantly-related organisms this metric is less useful. Viruses evolve very quickly, so AAI is more useful here. To get the mean percent identity we can use the "mean" function in R. 

> mean(blastout$V3)
 
I got 50.6 for this, so pretty low! This implies the two viruses are quite divergent (not particularly closely related). 
Note that this is the one-way AAI, since we are only looking at results using one virus as the query and one as the reference. Results will vary slightly if we do the reverse. 

Try doing the reverse and seeing how similar the results are (i.e., using PSSM2 as the query and PSSM3 as the db).
When you vary the e-value what happens to the one-way AAI? Does this make sense?
What about query coverage? How does increasing that change the one-way AAI?


