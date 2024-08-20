## Introduction to HMMs ##
<br>


First let's download this folder:

>git clone https://github.com/faylward/hmm_introduction

and then 

>cd hmm_introduction


You will find a FASTA file called polb.faa. This contains sequences of family B DNA polymerases that are encoded by several different kinds of herpesviruses. 

And now let's get some stats on how many sequences are in this file, and how long they are:

>seqkit fx2tab -inl polb.faa

Today, we want to make a Hidden Markov Model (HMM). For this we will need to first contruct a global alignment of different homologus proteins. There are many different multi-sequence alignment programs out there. Here we will use Clustal Omega. Because it is nice to visualize these things we will go to the main webpage and use the web interface:

If you have not installed clustal omega yet, you can do so with:

>conda install clustalo -c bioconda

Now to make the alignment we can use a simple command:

>clustalo -i polb.faa -o polb.aln

And now let's take a look at the file:

>more polb.aln

Note that the alignment file is still in FASTA format, though now there are "-" characters to signify gaps. FASTA format is generally used for alignments, though for visualization tools like MView are nicer since you can see all of the aligned regions more clearly. 

Let's see how long these sequences are:

>seqkit fx2tab -iln polb.aln

Note that all of the sequences are the same length now. This should always be the case for an alignment, since gap characters should effectively lengthen shorter sequences.

We can also visualize the alignment using some online tools. A simple one is appropriately called ALIGNMENTVIEWER: https://alignmentviewer.org/

Once you go to this site you can upload the alignment by clicking on the Upload button on the upper right hand side of the screen. If you are using an HPC for this tutorial, you will need to download the alignment onto your computer first.  

Browse through the alignment and note how some regions are highly conserved while others are highly variable. Which regions do you think will be most informative for classifying new sequences? 

Another alignment viewer is one offered by NCBI - you can access it here: https://www.ncbi.nlm.nih.gov/projects/msaviewer/


Now let's create an HMM from the multi-sequence alignment of PolB proteins. For this we can use the hmmbuild command in HMMER. 

>hmmbuild polb.hmm polb.aln

You should see the polb.hmm file now. Take a look with "more".  


Now we want to test this new polb.hmm out to see how well it predicts polb proteins.

Let's start with a well-characterized herpesvirus - the famous Epstein-Barr virus. We can download the genome of EBV with the following command:

> wget -O epstein_barr.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.fna.gz

and unzip it:

>gunzip epstein_barr.fna.gz

We downloaded a genome file, so to examine the encoded proteins we first need to predict protein sequences with Prodigal:

To compare a protein file to a HMM we can use the hmmsearch command in HMMER. 

> prodigal -i epstein_barr.fna -a epstein_barr.faa

Now we can search the proteins against our HMM using hmmsearch:

>hmmsearch polb.hmm epstein_barr.faa

Or, if we want a tabulated output and introduce an E-value threshold, 

>hmmsearch -E 1e-10 --tblout ebv_vs_polb.hmmout polb.hmm epstein_barr.faa

We can  browse the results with "more". Are there any good matches to our PolB HMM?

Now let's compare this annotation to one that we would do with BLASTP. 
First we need to make a BLASTP database from the polb.faa file. 

>makeblastdb -in polb.faa -dbtype prot

And then we can compare all of the EBV proteins to this BLAST database, using an e-value of 1e-10. 

>blastp -query epstein_barr.faa -db polb.faa -outfmt 6 -max_target_seqs 1 -evalue 1e-10

In the case of EBV, we get pretty good results with both BLASTP and hmmsearch, meaning that both tools recover the EBV PolB protein rather well. Let's try with another herpesvirus that is lesser known and is not closely related to EBV. Let's try a herpesvirus that infects abalone. 

First download the genome: 

>wget -O abalone_herpesvirus.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/900/375/GCF_000900375.1_ViralProj177933/GCF_000900375.1_ViralProj177933_genomic.fna.gz

Unzip it: 
>gunzip abalone_herpesvirus.fna.gz

And predict proteins: 

>prodigal -i abalone_herpesvirus.fna -a abalone_herpesvirus.faa

And run the hmmsearch: 

>hmmsearch --tblout abalone_vs_polb.hmmout polb.hmm abalone_herpesvirus.faa

Do we get a good hit to PolB here? 
Now let's cross-reference with the BLASTP results:

>blastp -query abalone_herpesvirus.faa -db polb.faa -outfmt 6 -max_target_seqs 5 -evalue 1e-10

Now the results here are less clear-cut. HMMs can recover the abalone herpesvirus PolB reasonably well, but BLASTP provides only a few hits with low percent identity. For the analysis of divergent protein families, it is typically best to use HMMs. 



