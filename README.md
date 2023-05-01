# entropy
## Scripts for entropy calculation and filtering

# shannons-filter.py

Theo Allnutt 2021

r3.7 Calculates Shannon's information index for sample of reads, only outputs reads over threshold.
Usage:

shannons-filter.py "reads/*" filtered/ 0.8 12 paired 20000000

filtered/=output folder, 0.8= threshold, 12=threads, paired=paired reads, 20000000 number of reads to sample

Accepts fasta, fastq and .gz
