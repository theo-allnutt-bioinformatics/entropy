# Entropy
## Scripts for entropy calculation and filtering

## shannons-filter.py

Theo Allnutt 2021  
r3.7 Calculates Shannon's information index for sample of reads, only outputs reads over threshold.  

Usage:  

	shannons-filter.py "<read folder>/*" <filtered folder>/ <threshold> <threads> <paired/single> <max reads to sample>  
	shannons-filter.py "reads/*" filtered/ 0.8 12 paired 20000000  

filtered/=output folder, 0.8= threshold, 12=threads, paired=paired reads, 20000000 number of reads to sample  

Accepts fasta, fastq and .gz


## kmercount-shannons.py

r10.2 outputs Shannon's and kmer count entropy measures for DNA sequence file or files.  

Usage:  

	kmercount-shannons.py "<input-folder>/*" <output-folder>/ <kmerlen> <threads>  
	kmercount-shannons.py "mega353.fasta" ./ 8 1  
	kmercount-shannons.py "reads/*.fasta" kmercounts/ 8 12  
	kmercount-shannons.py "reads/*.fastq.gz" kmercounts/ 4 12  
