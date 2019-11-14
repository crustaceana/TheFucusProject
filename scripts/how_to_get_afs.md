# What  
The script `get_afs.R` can be executed in order to perform two main tasks:

1. Determine ancestral vs. derived state of alleles.  
The outgroups are two sexual populations, "KRI_sex" and "TJA_sex", and the allele is ancestral when the variant is monomorphic in both the outgroup populations. Obviously, the other allele is the derived one. We are exclusively analysing biallelic variants.

1. Plot derived allele frequency of one population against a second population.

# How  
To run `get_afs.R` the user must be in the same directory where the `data/` directory is (this directory is not on GitHub). The present [repository](https://github.com/crustaceana/TheFucusProject.git) must also be in the same place (i.e., if you run `ls` on your terminal you should see at least something like: `data/` `TheFucusProject/`).  
Copy and paste the command below to launch the script:  
`Rscript TheFucusProject/scripts/get_afs.R -C data/YOURDATA.csv -E .svg -O "KRI_sex TJA_sex"`

After pressing return, some R packages will be loaded and those that do not already exist will be first installed and then loaded. If no has occurred, the next step of the analysis is to calculate the allele frequencies and your screen will be filled with these type of lines:

```
...
Calculating frequency of scaffold13248_99885 ...
Calculating frequency of scaffold13249_30878 ...
Calculating frequency of scaffold13249_40170 ...
Calculating frequency of scaffold13262_28955 ...
Calculating frequency of scaffold13262_35832 ...
Calculating frequency of scaffold13263_17446 ...
...
```

Not very useful but it is just a check for the accuracy of the computation. Much more useful are the lines towards the end of the analysis:

```
The number of polarised variants is 543
There are 2 variants that have non-zero frequency in all populations.
There are 182 variants that have non-zero frequency in the sexual populations.
The correlation between LET_sex and STO_sex frequencies is 0.3383625
There are 2 variants that have non-zero frequency in all the asexual populations.
The variants with non-frequency that are shared between sexual and clonal populations are:
[1] "AF_scaffold10603_55680" "AF_scaffold12082_49062"
```

Other outputs will appear on the screen but again these are mainly check points.

PS: If you cannot guess what the flags of the Rscript command mean there is a help flag which will return all the details:  
`Rscript TheFucusProject/scripts/get_afs.R -h`
