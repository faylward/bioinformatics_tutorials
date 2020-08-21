
Download the VOG database

> mkdir vog
> cd vog
> wget http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz
> cd ..
> cat vog/*.hmm > vogdb.hmm

Run the VOG annotation

> python vog_annotation.py proteins/ vog_table.txt

