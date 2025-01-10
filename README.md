# General
This pipeline is for bacterial transcriptomes. It takes raw reads of ONT direct RNA sequencing (or stranded long reads) along with the genome reference and map the full features of the trasncriptome. This included trasncription start and termination sites, operon maps, reconstructed trasncriptome with focus on anti-sense and intergenic RNA.

Output files header:
TSSs:  Gene name, TSS, number of reads supporting this TSS (number of long reads that start around this site)

Operons: first gene-last gene(or operon name), all genes in the operon, number or reads supporting this operon (~ count of reads corresponding to this operon), relative expression (for example if operon A-D count is 30 reads and gene A is 50 then relative expression is 3/5 = 0.6), relative expression for all genes in the operon in order.

Final transcripts: 
Chromosome, start site, end site, name, strand, count, type, overlapping genes (from reference annotation)
![image](https://github.com/user-attachments/assets/dd0ad70e-19b1-4313-ae82-093d075662d0)
