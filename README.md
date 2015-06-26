# ono
Order and orient

This tool tries to order and orient the contigs of a genome using a reference.
It is based on the libraries of : https://github.com/jasperlinthorst/reveal


	python assembly_finishing.py --minmum minimal_match_size -n number_of_N_to_insert_between_contigs reference.fasta contigs.fasta
optional : 
	-discard don't append the contigs that couldn't be ordered and orientated to the output
	-p prune all sequences of N to a single N (overrides the -n option), even those present in the input contigs