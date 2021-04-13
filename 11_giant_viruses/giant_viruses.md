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

let's download the genome for Acanthamoeba polyphaga Mimivirus, the canonical "giant virus":
> wget -O mimi.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/888/735/GCF_000888735.1_ViralProj60053/GCF_000888735.1_ViralProj60053_genomic.fna.gz

and unzip it

> gunzip mimi.fna.gz

and then run VR on it:

python viralrecall.py -i mimi.fna -p mimi_out -f

