# What  
The script ´get_afs.R´ can be executed in order to perform two main tasks:

1. Determine ancestral vs. derived state of alleles.  
The outgroups are two sexual populations, "KRI_sex" and "TJA_sex", and the allele is ancestral when the variant is monomorphic in both the outgroup populations. Obviously, the other allele is the derived one. We are exclusively analysing biallelic variants.

1. Plot derived allele frequency of one population against a second population.

# How  
To run ´get_afs.R´ the user must be in the same directory where the ´data/´ directory is (this directory is not on GitHub). Also, the present [repository](https://github.com/crustaceana/TheFucusProject.git) must also be in the same place (i.e., if you run ´ls´ on your terminal you should see at least something like: ´data/´ TheFucusProject/´).  
Copy and paste the command below to launch the script:  
´Rscript TheFucusProject/scripts/get_afs.R -C data/YOURDATA.csv -T YOUROUTPUT.txt -O "KRI_sex TJA_sex"´

If you cannot guess what the flags mean there is a help flag which will return all the details:  
´Rscript TheFucusProject/scripts/get_afs.R -h´
