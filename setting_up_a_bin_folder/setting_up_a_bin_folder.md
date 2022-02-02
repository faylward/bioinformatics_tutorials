
#### Setting up a bin folder #####

Why set up a bin folder on your unix machine (or home folder on an HPC) ?

If we download a binary we would like to be able to execute it from the command line regardless of what directory we happen to be in. 
For this it is convenient to set up a bin folder where we can put all of the binaries that we want to download, and then configure the command line such that it looks there every time we execute a command. 

First, we need to create the bin folder. Navigate to your home directory and then run:

>mkdir bin

Then we need to find a binary to download and put in there. Let's start with Seqkit, a handy tool for manipulating FASTA files. 
Let's take a look at the download page: https://bioinf.shenwei.me/seqkit/download/

From here we can find the URL for the appropriate download- on my Unix system it is the linux_amd64 version. We can navigate to the bin folder and download this. 

> cd bin

and then:

>wget https://github.com/shenwei356/seqkit/releases/download/v2.1.0/seqkit_linux_amd64.tar.gz

And then we need to unzip it, which in this case would require a tar command. 

> tar xvfz seqkit_linux_amd64.tar.gz

And then we can check to make sure it is there. 

> ls

Now we should be able to run seqkit in the command line as long as we are in the bin folder. But if we configure a the .bashrc file in our home directory, we will be able to run binaries in our bin folder from any directory. 
The .bashrc file in your home directory contains code that is run every time the command line is initiated.  



and

>cd transcriptomics_tutorial
>
and then

>ls
