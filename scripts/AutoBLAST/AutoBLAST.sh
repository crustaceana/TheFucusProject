#! /bin/bash

# AutoBLAST
# Matt Pinder, 2020

# Given a list of SNPs and a reference genome, BLAST the 10kb surrounding the SNP (5kb up- and downstream)

# Requires BLAST and Python2

# Arguments:
# -s = List of SNPs, one per line, in format contig_position
# -r = Reference genome
# -t = Number of threads for BLAST (optional)


# If number of threads not specified, then use single thread
THREADS=1

while getopts :s:r:t: opt; do
	case $opt in
		s)
			SNPLIST=$OPTARG
			echo "List of SNPs: $OPTARG"
		;;
		r)
			REF=$OPTARG
			echo "Reference: $OPTARG"
		;;
		t)
			THREADS=$OPTARG
			echo "Number of threads: $OPTARG"
		;;
	esac
done

####################

while read i; do
	echo "Current SNP: $i"

	# Determine the 'contig' and 'position' parameters for each SNP
	CONTIG=$(echo ${i} | cut -f1 -d'_')
	POSITION=$(echo ${i} | cut -f2 -d'_')

	# Determine the upper and lower boundaries to cut between (5kb up- and downstream of the SNP)
	# Note that cut can handle too-high positions, but not negatives, hence the if-statement
	if [ $POSITION -le 5000 ]; then
		START=1
	else
		let "START = $POSITION - 5000"
	fi
	let "END = $POSITION + 5000"

	# Use fp.py (by Mats Töpel; https://github.com/mtop/ngs) and Bash to extract the desired region from the reference genome
	QUERY=$(fp.py --seq $CONTIG $REF | tail -n1 | cut -c${START}-${END})

	# Run BLASTn on the query (feel free to adjust this line to suit your own BLASTing needs)
	echo -e ">${CONTIG}:${START}-${END}\n${QUERY}" | blastn -db /db/nt -num_threads $THREADS -out ${i}_vs_nt.BLASTN.txt

done < "$SNPLIST"

echo "Done!"
