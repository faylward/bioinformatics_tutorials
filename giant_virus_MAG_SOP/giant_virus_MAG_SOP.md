## Example workflow for binning and analyzing giant virus MAGs!

Several people have asked me to write up a short workflow on how we produce giant virus MAGs in our lab. 
What follows here is a fairly bare-bones workflow, and bells-and-whistles can be added to suite the particular needs of different datasets. 
The steps below outline methods we use in the lab for:

-Producing giant virus MAGs

-Quality-checking MAGs

-Performing functional annotations on proteins encoded in the MAGs

-Placing the MAGs into a reference phylogeny for classification. 

This is based on a workflow we developed in the lab for a previous study (Moniruzzaman et al., Nature Communications, 2020). 

Here I am using the term "giant virus" to refer to all dsDNA viruses in the Nucleocytoviricota - a viral phylum that includes well-studied viruses such as 
poxoviruses and African swine fever virus as well as many more-recently discovered lineages such as Pandoraviruses and Mimiviruses (among many others). 
The nomenclature is really confusing and many quasi-contradictory taxonomic names have been proposed over the years - for more details on how we are trying to 
standardize things please see our preprint on giant virus taxonomy: 

It turns out giant viruses are quite ubiquitous in the biosphere, and they turn up in a wide variety of metagenomes. So I anticipate that
workflows similar to this one could be applied on a wide variety of datasets, provided that assembly quality was high enough. As a good rule-of-thumb, 
if a metagenome is high enough quality to recover bacterial or archaeal MAGs it should be good enough to do the same for giant viruses. 

But first of all, why would we bother generating viral bins in the first place? Most metagenome-based viral diversity studies don't do this - 
instead they rely on direct analysis of viral contigs that have been retrieved from an assembled metagenome. This contig-centric approach might work alright
for smaller viral genomes (i.e., some bacteriophage), but giant viruses have large and complex genomes - typically 200-1,200 kbp - and they will usually not
assemble down into nice individual contigs. So by using a contig-centric approach we would be missing lots of valuable information since multiple contigs
belong to the same virus. So in my opinion it is best to analyze giant viruses in a similar manner to how we routinely analyze bacteria and archaea - by 
generating MAGs using binning-based approaches and subsequent QC. 

For binning purposes I am just going to use MetaBat2 here. I like MetaBat2 because it is scalable and fairly easy to use, but please be advised that 
there are many different approaches for contig binning out there, and there may be a better approach for your data. 
You may wish to use multiple binning approaches and see which one works best. For example if you had multiple metagenomes from 
similar environments you may wish to incorporate contig co-abundance information into your analysis. I'm not going to go into all the different binning
strategies out there because it would lead to a very long digression. The main principle to keep in mind is that coverage and nucleotide composition (i.e., 
tetranucleotide frequencies) can be used to bin giant virus genomes just like they can be used to make bacterial/archaeal MAGs. Just make sure the 
approach you use does not use marker genes, because giant viruses will not have the same marker genes we would expect to find in bacteria/archaea. 

download contigs
download raw reads
map reads to get coverage information
use MetaBat2 to get bins
run ViralRecall to identify NCLDV bins
Place NCLDV bins into the reference GVDB phylogeny
Analyze annotation results from ViralRecall




first we

>mkdir bins

>metabat2 -i sample_data/final.contigs.fa -a sample_data/HOT229_0025m.coverm -o bins/testbin -s 100000 -t 4 -m 10000

This is just a simple example with default parameters. I am being a bit cautious here and only considering contigs >10 kbp in length. 
If you want you can play around with the parameters a bit- in the past I have made things a bit more stringent with things like "--minS 75 --maxEdges 75"", but in my experience this rarely changes much (maybe a single small contig will be included in with the default paramters but excluded with more stringent ones).

Now that we have bins we can see which ones belong to giant viruses. For this #I'll use ViralRecall, which I developed. 

>git clone https://github.com/faylward/viralrecall

and move the data into this directory for simplicity:

>mv bins viralrecall
>cd viralrecall

We need to download and unpack the appropriate reference HMM databases in order to run ViralRecall. 

>wget -O hmm.tar.gz https://zenodo.org/record/4762520/files/hmm.tar.gz?download=1

and then unpack:
>tar -xvzf hmm.tar.gz

Now we can run ViralRecall on the bins. If you have never used ViralRecall you will want to make sure you have HMMER3 and Prodigal installed in your PATH, and that you have various python modules installed: Biopython, matplotlib, pandas, numpy, and defaultdict (if you are not sure if you have these installed you can just try running it and Python will tell you which ones you need to install- you can usually do this with Pip). I won't go over all the details here, but more is explained on the ViralRecall github page (https://github.com/faylward/viralrecall).

>python viralrecall.py -i bins -p bins_out -b -c -f -t 4

The -b flag tells the program to run iteratively over all files in a directory. I'm also going to use the -c flag so that the results are provided on a per-contig basis (not a rolling average across a contig, which would be more useful for identifying endogenous viruses), and I'll use the -f flag to produce plots for each of the bins. This command took a little over an hour to run on my system.

ViralRecall will not only tell us which bins belong to giant viruses, it will also predict proteins using Prodigal and perform a full protein-level annotation using the Pfam and Giant Virus Orthologous Groups (GVOG) databases. But because ViralRecall does so much, it is also a bit slow, so if you have over a few hundred bins you may wish to set up multiple jobs in parallel. If you have lots of MAGs there are also a few alternative options you may consider:

Alternative A) If you were interested in bacterial/archaeal bins as well you could run your bins through CheckM/GTDB-tk and separate all the bins that have high completeness or a good classification there- that way you wouldn't need to run ViralRecall on bins that you already know belong to bacteria or archaea. So however you organize this will depend on your overall goals and the size of your dataset. 

Alternative B) You can run ViralRecall with the "-db marker" option, which will only run a HMMER search against a small set of giant virus marker genes and is therefore much faster. You can look at the results in the *.summary.tsv files to see what markers are found on which contig. If you do this just beware that some markers are also found in bacteria/archaea (the RNA polymerase subunits, for example), so if you will want to look for the markers that are more specific to giant viruses - the major capsid protein (MCP), VLTF3 transcriptional factor, and A32 packaging ATPase. If you have bins where some of these markers are present you can opt for a full ViralRecall run, but if no markers are present you can save some time and avoid analyzing the bins further. 


Anyway if we run ViralRecall on all the bins here with the command above we get a folder called bins_out with all the results. In this folder is a file called batch_summary.txt that gives us an idea of what the results look like. We are looking for bins where all or almost all of the contigs have positive scores (more similar to giant viruses) - based on this criteria we have four promising bins here- testbin 5, 6, 10, and 13. To look at these in more detail we can look at the *.summary files and .pdf figures that were generated. 



