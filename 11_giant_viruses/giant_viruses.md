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

