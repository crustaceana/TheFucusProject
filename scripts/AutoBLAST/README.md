# Script for extracting the region around a SNP of interest and BLASTing it,  
  and the SGE script required to run it on the Albiorix cluster

## Usage:
`./AutoBLAST.sh -s [SNP list] -r [Reference genome] -b [BLAST database] -t [Number of threads]`

## Notes
* Requires `fp.py` (see https://github.com/mtop/ngs), Python2 and BLAST in your path
* The SNP list must be in the format `contig_position`
* If not being run on the Albiorix cluster, or if different BLAST settings are used,  
  the BLAST command will need to be edited
