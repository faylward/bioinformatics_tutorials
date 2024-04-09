

### Intro to metagenomic classification ###

First we need to install kraken2

> conda install kraken2 -c bioconda

Then we need to download the appropriate database from the kraken2 website

> wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz

And unpack it:

> tar xvfz k2_standard_08gb_20240112.tar.gz

Then we need to download some data to start processing, using our favorite fastq-dump command. I happend to find this metagenome before - it is from seawater. We will only download 20,000 sequences to start. 

> fastq-dump -X 20000 --split-3 SRR5131972

And now we can run kraken2 to classify the reads in this metagenome! For best results, make sure to say "release the kraken" when running this tool. 

> kraken2 --use-mpa-style --report SRR5131972.report.txt --output SRR5131972.classification.txt --db k2_standard_08gb SRR5131972_1.fastq SRR5131972_2.fastq

We can see the classification for each read in the *classification.txt file, and we can see the summarized results in the *report.txt file. 

Now, it would be nice if we could make a Krona diagram that summarizes this data for us. For that we will need to install a tool called recentrifuge
For that we can use pip.

> pip install recentrifuge

And hopefully it works - the command is "rcf"

> rcf -h

This tool needs NCBI taxonomy files to run properly, so we can download those files with the following command:

> retaxdump 

You will want to move those files into a directory called "taxdump" and then unzip it:

> mkdir taxdump
> mv taxdmp.zip taxdump
> unzip taxdump/taxdmp.zip

And now we can run recentrifuge

> rcf -k SRR5131972.classification.txt --sequential

Now let's look at a metatranscriptome from the same sample!

> fastq-dump -X 20000 --split-3 SRR5131969
> kraken2 --use-mpa-style --report SRR5131969.report.txt --output SRR5131969.classification.txt --db k2_standard_08gb SRR5131969_1.fastq SRR5131969_2.fastq
> rcf -k SRR5131969.classification.txt --sequential
