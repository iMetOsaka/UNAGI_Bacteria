#Unagi for bacteria
import sys, os, argparse, subprocess


from log import logger
import confb
log=logger()
config=confb.getconf()
version="v1.0.0"


def main(argv):
	global log
	global config
	global version
	## Arguments definition ##
	parser = argparse.ArgumentParser("UNAGI pipeline.")
	parser.add_argument('-v', '--version', action='version', version="%(prog)s ("+version+")")
	parser.add_argument('-i', '--input', help='Input file. This should be the path to a fastq file.', required=True)
	parser.add_argument('-o', '--output', help='Output path. This should be the path to output the result files to.', required=True)
	parser.add_argument('-a', '--annotation', help='Annotation file. It should be in gff', required=False)
	parser.add_argument('-g', '--genome', help='Genome File. It should be in fasta', required=True)	
	parser.add_argument('-V', '--verbose', help='Verbose. Displays additional information while running.', action='store_true', required=False)
	parser.add_argument('-S', '--silent', help='Runs silently (no console output), cancels verbose.', action='store_true', required=False)
	args = parser.parse_args(argv)

	#Config file check
	if config == None:
		log.tell("The configuration file is not ini formatted")
		return

	## Arguments check ##
	#Required
	#Input
	inputFile = args.input
	if not os.path.isfile(inputFile):
		log.tell("The input file %s doesn't exist"%(args.input))
		return
	#Output
	outputPath = args.output
	if not os.path.isdir(outputPath):
		log.tell("The output directory %s doesn't exist, creating it"%(outputPath))
		os.mkdir(outputPath)

	transitionnalOutputPath=os.path.join(outputPath,"transitionnal")
	if not os.path.isdir(transitionnalOutputPath):
		os.mkdir(transitionnalOutputPath)


	#Genome
	geneFile = args.annotation
	if not os.path.isfile(geneFile):
		log.tell("The annotation file %s doesn't exist"%(args.annotation))
		return
	genomeFile = args.genome
	if not os.path.isfile(genomeFile):
		log.tell("The genom file %s doesn't exist"%(args.genome))
		return


	#More output should happen if the script runson verbose mode
	log.verbose = args.verbose
	#No output should happen if the script runs silently
	log.silent = args.silent
	
	#files names
	samFile = os.path.join(transitionnalOutputPath,config["raw_mapped_sam_file"])
	bamFile= os.path.join(transitionnalOutputPath,config["raw_mapped_bam_file"])
	chrFile = os.path.join(transitionnalOutputPath,config["chr_file"])
	bamsortedFile = os.path.join(transitionnalOutputPath,config["raw_sorted_bam_file"])
	bedFile = os.path.join(transitionnalOutputPath,config["raw_sorted_bed_file"])
	awarFile = os.path.join(transitionnalOutputPath,config["awar_tss_file"])
	intermediateFile,unawarFile = os.path.join(transitionnalOutputPath,config["clustering_TTS_raw"]),os.path.join(transitionnalOutputPath,config["unawar_tss_file"])
	rawTssFile,finalTssFile = os.path.join(transitionnalOutputPath,config["raw_list_TSSs_file"]),os.path.join(transitionnalOutputPath,config["final_list_TSSs_file"])
	poscovFile,negcovFile = os.path.join(transitionnalOutputPath,config["positive_coverage_file"]),os.path.join(transitionnalOutputPath,config["negative_coverage_file"])
	posbedFile,negbedFile = os.path.join(transitionnalOutputPath,config["positive_bed_file"]),os.path.join(transitionnalOutputPath,config["negative_bed_file"])
	geneBedFile,intersectFile=os.path.join(transitionnalOutputPath,config["genes_bed_file"]),os.path.join(transitionnalOutputPath,config["intersect_file"])
	TTSFile=os.path.join(outputPath,config["list_TTSs_file"])
	intersectForOperonFile =os.path.join(transitionnalOutputPath,config["intersectforoperon_file"])
	operonsFile =os.path.join(outputPath,config["operons_file"])
	rawTrasncriptsFile =os.path.join(transitionnalOutputPath,config["raw_trasncripts_file"])
	trIntesectFile,rawBedTrFile,rawBedTrFilesorted,FinalTranscriptsFile=os.path.join(transitionnalOutputPath,config["transcripts_intersect_file"]),os.path.join(transitionnalOutputPath,config["trasncripts_bed_file"]),os.path.join(transitionnalOutputPath,config["trasncripts_sortred_bed_file"]),os.path.join(outputPath,config["Final_trasncripts_file"])
	#Map the reads to the genome
	log.tell("Mapping the reads to the genome")
	minimap(inputFile, genomeFile,samFile )

	#Get a bam file from the results and sort it
	log.tell("Sorting the mapped reads")
	samToBam(samFile,bamFile,chrFile)
	sortBam(bamFile,bamsortedFile)
	
	#converting to bedfile
	log.tell("Converting to bed file")
	bamToBed(bamsortedFile,bedFile)

	#Identifying TSSs in an awar approach
	log.tell("Identigying TSSs")
	intersect(bedFile,geneFile,geneBedFile,intersectFile,intersectForOperonFile)
	AwaridentifyTSS(intersectFile,awarFile)
	#Identifying TSSs in an unawar approach
	log.tell("Identigying TSSs in anawar approach")
	seperateBedandGenomcov(bedFile,posbedFile,negbedFile,poscovFile,negcovFile,chrFile)
	UnAwarindentifyTSS(chrFile,bedFile,intermediateFile,unawarFile,poscovFile,negcovFile,rawTrasncriptsFile)
	#combining the two lists
	log.tell("Combining both TSSs lists")
	combineTSSs(unawarFile, awarFile, rawTssFile,finalTssFile,chrFile)
	#Identifying TTSs
	log.tell("Identifying TTSs")
	indentifyTTS(intersectFile, TTSFile)
	#identifying operon
	log.tell("Identifying Operons")
	identifyOperons(intersectForOperonFile,operonsFile)
	#constructing trasncriptome
	reconstruct(geneBedFile,rawTrasncriptsFile,trIntesectFile,rawBedTrFile,rawBedTrFilesorted,FinalTranscriptsFile)


def intersect(Bedfile,geneFile,geneBedFile,intersectFile,intersectForOperonFile):
	global config
	global log
	log.write("Converting GFF to Bed")
	type = config["gene_name"]
	with open(geneFile,'r') as gff,open(geneBedFile,'w') as new:			
		for line in gff: #store genes coordinates first
			read = line.strip().split('\t')
			#first checking if it is a gene and extracting its name from a gff file
			if line[0] != '#' and read[2] == "gene":
				start, end = int(read[3]), int(read[4])
				strand = read[6]
				chromosome = read[0]
				name = ""
				last = read[-1].split(';')
				for i in last:
					if i[0:len(type)] == type:
						name = i[len(type)+1:]
						break
				if name != "":
					new.write('%s\t%i\t%i\t%s\t60\t%s\n'%(chromosome,start,end,name,strand))				
	command=config["bedtools_path"]+" intersect -wo -a " +geneBedFile +" -b " + Bedfile + " > "+ intersectFile
	log.write("Running bedtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	command=config["bedtools_path"]+" intersect -wo -a " +Bedfile +" -b " +geneBedFile  + " > "+ intersectForOperonFile
	log.write("Running bedtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running bedtools")    
def bamToBed(sourceFile, outputFile):
	global config
	global log
	command=config["bedtools_path"]+" "+config["bamtobed_options"]+" "+sourceFile+" > "+outputFile
	log.write("Running bedtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running bedtools")
def AwaridentifyTSS(intersectFile,outputFile):
	global config
	global log
	threshold = int(config["min_TSS_threshold"])
	cov_threshold = float(config["min_coverage_for_TSS"])
	
	with open(intersectFile, 'r') as gen, open(outputFile, 'w') as new:
		new.write('Chromosome\tTSS\tCount\tGene\tstrand\tPredictedTTS\n')
		first = gen.readline()
		read = first.strip().split()
		gene = read[3]
		sign = read[5]
		start = int(read[1])
		end = int(read[2])
		chromosome = read[0]
		all_overlap = [(int(read[7]),int(read[8]))]
		for line in gen:
			read = line.strip().split()
			if read[3] == gene:
				if read[11] == sign:
					all_overlap.append((int(read[7]),int(read[8])))
			else:
				if sign =='+':
					after = list()
					all_overlap.sort()
					first_start = all_overlap[0][0]
					index = -1
					if first_start < start:#smallest start is before start codon
						while first_start < start and index < len(all_overlap)-1:
							index +=1
							first_start = all_overlap[index][0]
						if index < len(all_overlap):
							after = all_overlap[0:index+1]
						else:
							after = all_overlap
						#now resolve if there are multiple peaks
						current = list()
						ends = list()
						cov = len(after) # the coverage at the start codon is the number of overlapping reads
						after.sort() #sort again according to the start site
						prev = after[0][0]
						for end in after: #end is actually the read not just an end
							if end[0] - prev < threshold: #if the distance between the read start site and the next one is less than the threshold, take in the same cluster
								current.append(end[0])
								ends.append(end[1])
								prev = end[0]
							else:
								if len(current) > (cov_threshold * cov):  # if the read is far from than the next one, finish the current cluster and check it if its count is important
									ends.sort() # optional step to calculate the end of these start sites
									tss_end = ends[-1] # the end is the most extrem end
									new.write('%s\t%i\t%i\t%s\t+\t%i\n' % (chromosome, current[0] + 1, len(current), gene, tss_end))
								current = [end[0]] #start a new cluster with the current read
								ends = [end[1]]
								prev = end[0]
						if len(current) > (cov_threshold * cov) and len(current) > 2: #to check the last cluster
							ends.sort()
							tss_end = ends[-1]
							new.write('%s\t%i\t%i\t%s\t+\t%i\n' % (chromosome, current[0] + 1, len(current), gene, tss_end))
					else:
						new.write('%s\t%i\t%i\t%s\t+\t%i\n' % (chromosome, all_overlap[0][0], 0, gene, all_overlap[-1][1]))
					gene = read[3]
					sign = read[5]
					start = int(read[1])
					end = int(read[2])
					chromosome = read[0]
					all_overlap = [(int(read[7]),int(read[8]))]
				else:
					
					all_overlap.sort(key=lambda all_overlap: all_overlap[1])
					first_start = all_overlap[-1][1]
					if first_start > end:
						before = list()
						for site in reversed(all_overlap):
							if site[1] < end:
								break
							else:
								before.append((site[1],site[0]))
						before.sort()
						prev = before[0][0]
						cov = len(before)
						current = list()
						ends = list()
						for i in before:
							if (i[0] - prev) < threshold:
								current.append(i[0])
								prev = i[0]
								ends.append(i[1])
							else:
								if len(current) > (cov_threshold * cov):
									ends.sort()
									tss_end = ends[0]
									new.write('%s\t%i\t%i\t%s\t-\t%i\n' % (chromosome, current[-1], len(current), gene, tss_end))
								current = [i[0]]
								prev = i[0]
						if len(current) > (cov_threshold * cov) and len(current) > 2:
							ends.sort()
							tss_end = ends[0]
							new.write('%s\t%i\t%i\t%s\t-\t%i\n' % (chromosome, current[-1], len(current), gene, tss_end))
					gene = read[3]
					sign = read[5]
					start = int(read[1])
					end = int(read[2])
					chromosome = read[0]
					all_overlap = [(int(read[7]),int(read[8]))]
def UnAwarindentifyTSS(chrFile,bedFile,intermediateFile,outputFile,pos_cov,neg_cov,rawTrasncriptsFile):
	global config
	global log
	threshold = int(config["min_TSS_threshold"])
	chromosomes = set()
	with open(chrFile,'r') as chr:
		for line in chr:
			read = line.split()
			chromosomes.add(read[0])		
	#input file is a bed file generated from sorted bam alignment file using bedtools bamtobed
	with open(bedFile, 'r') as old, open(intermediateFile, 'w') as new, open(rawTrasncriptsFile,'w') as tr:
		pos_reads =dict() #seperate positive from negative reads
		neg_reads=dict()
		for chromosome in chromosomes:				
			pos_reads[chromosome] = list() #seperate positive from negative reads
			neg_reads[chromosome] = list()
		for line in old:
			read = line.strip().split()
			if read[-1] == '+':
				pos_reads[read[0]].append((int(read[1]),int(read[2])))
			else:
				neg_reads[read[0]].append((int(read[2]),int(read[1])))
		for chromosome in pos_reads:
			chromo = pos_reads[chromosome]
			chromo.sort()
			prev = chromo[0][0]
			current = [chromo[0][0]]
			ends = [chromo[0][1]]
			for site in chromo[1:]:
				if site[0] - prev < threshold:
					current.append(site[0])
					ends.append(site[1])
					prev = site[0]
				else:
					if len(current) > 3: #if the cluster is more than three reads considere it for further analysis
						ends.sort()
						new.write('%s\t%i\t%i\t%i\t+\t%i\n' % (chromosome, current[0]+1, len(current), current[-1] - current[0],ends[-2]))
						clustered_ends = clustering(ends,50)
						for end in clustered_ends:
							tr.write('%s\t%i\t%i\t%i\t+\n' % (chromosome, current[0]+1, end[0],end[1]))
					prev = site[0]
					current = [site[0]]
					ends = [site[1]]

		for chromosome in neg_reads:
			chromo = neg_reads[chromosome]
			chromo.sort()
			prev = chromo[0][0]
			current = [chromo[0][0]]
			ends = [chromo[0][1]]
			for site in chromo[1:]:
				if site[0] - prev < threshold:
					current.append(site[0])
					ends.append(site[1])
					prev = site[0]
				else:
					if len(current) > 3:
						ends.sort()
						new.write('%s\t%i\t%i\t%i\t-\t%i\n' % (chromosome, current[-1], len(current), current[-1] - current[0],ends[1]))
						clustered_ends = clustering(ends,50)
						for end in clustered_ends:
							tr.write('%s\t%i\t%i\t%i\t-\n' % (chromosome, current[-1],end[0]+1,end[1]))					
					prev = site[0]
					current = [site[0]]
					ends = [site[1]]

	# filter according to coverage, it was calculated with bedtools genomecoverage
	with open(intermediateFile, 'r') as old, open(outputFile, 'w') as new, open(pos_cov, 'r') as pos, open(neg_cov, 'r') as neg:
		coverage = {'-':dict(),'+':dict()}
		for chromosome in chromosomes:
			coverage['+'][chromosome] = dict()
			coverage['-'][chromosome] = dict()
		for line in pos:
			read = line.strip().split()
			coverage['+'][read[0]][int(read[1])+1] = int(read[2])
		for line in neg:
			read = line.strip().split()
			coverage['-'][read[0]][int(read[1])+1] = int(read[2])
		for line in old:
			read = line.strip().split()
			if int(read[2]) > coverage[read[-2]][read[0]][int(read[1])]:  # if cluster meet the coverage criteria leave it, else don't print it
				new.write(line)
def combineTSSs(unawar,awar,rawOutput,outputFile,chrFile):
	chromosomes = set()
	with open(chrFile,'r') as ch:
		for line in ch:
			read = line.strip().split()
			chromosomes.add(read[0])
	with open(awar,'r') as aw, open(unawar,'r') as unaw, open(rawOutput,'w') as new:
		all_sites = dict()
		for chromosome in chromosomes:
			all_sites[chromosome] = list()
		aw.readline() #to pass the header
		for line in aw:
			read = line.strip().split()		
			all_sites[read[0]].append((int(read[1]),line))
		for line in unaw:
			read = line.strip().split()			
			all_sites[read[0]].append((int(read[1]),line))		
		for chr in chromosomes:
			chromosome = list(all_sites[chr])
			chromosome.sort()
			for site in chromosome:
				new.write(site[1])
	with open(rawOutput,'r') as old, open(outputFile,'w') as new:
		common = 0
		counts = dict()
		annotated = dict()
		ends = dict()
		for line in old:
			read = line.strip().split('\t')
			if not read[3].isdigit(): #our gene name all start with gene, if it is not the case you can use if not read[3].isdigit()
				name = read[0] + ';' + read[1] + ';' + read[4]
				annotated[name] = read[3]
			else:
				name = read[0] + ';' + read[1] + ';' + read[4]
			if name in counts:
				common +=1
				if counts[name] < int(read[2]):
					counts[name] = int(read[2])
					ends[name] = int(read[5])
			else:
				counts[name] = int(read[2])
				ends[name] = int(read[5])
		for site in counts:
			end = ends[site]
			if site in annotated:
				name = annotated[site]
			else:
				name = 'novel'
			read = site.split(';')
			new.write('%s\t%i\t%i\t%s\t%s\t%i\n'%(read[0],int(read[1]),counts[site],read[2],name,end))
		print('number of duplicated TSSs is %i'%(common))
def indentifyTTS(intersectFile, outputFile):
	global config
	global log
	#identification TTS by intersection
	threshold = int(config["min_TSS_threshold"])
	tss = list()

	#input file is obtained by bedtools intesect -a genes_file.bed -b reads_file.bed -wo
	with open(intersectFile,'r') as bed, open(outputFile,'w') as new:
		prev_gene = ''
		current = list()
		for line in bed:
			read = line.strip().split()
			if read[3] == prev_gene:
				if read[5] == read[11] and (min(int(read[2]),int(read[8]))-max(int(read[1]),int(read[7])))/(int(read[2])-int(read[1])) > 0.9:
					if read[11] == '+' and int(read[8]) > int(read[2]):
						current.append((int(read[8])))
					elif read[11] == '-' and int(read[1]) > int(read[7]):
						current.append((int(read[7])))
			else:
				if len(current) > 0 and prev_gene !='':
					current.sort()
					for site in clustering(current,threshold):
						new.write('%s\t%i\t%i\n'%(prev_gene,site[0],site[1]))
				prev_gene = read[3]
				current = []
		if len(current) > 0:
			current.sort()
			for site in clustering(current,threshold):
				new.write('%s\t%i\t%i\n' % (prev_gene, site[0], site[1]))
def identifyOperons(intersectFile,outputFile):
	global config
	global log
	coverage_threshold = float(config["operon_coverage_threshold"])
	quality_threshold = int(config["operon_readQuality_threshold"])
	if coverage_threshold > 1:
		print('Error, threshold should be less than 1')
	with open(intersectFile,'r') as old, open(outputFile,'w') as new:
		count_operon = dict()
		operon_name = dict()
		FirstLine = old.readline() #first line		
		read = FirstLine.strip().split('\t') 
		f_read = read
		prev = read[3]
		operon = [read[9]]
		sign = read[5]
		for line in old:
			read = line.strip().split()
			if read[3] == prev and read[5] == read[11] and int(read[4]) == quality_threshold and int(read[12]) > (
					coverage_threshold * (int(read[8]) - int(read[7]))):
				operon.append(read[9])
				#sign = read[5]
			else:
				if len(operon) > 1:
					#first read is not checked for sign
					if f_read[5] == f_read[11] and int(f_read[4]) == 60 and int(f_read[12]) > (
							coverage_threshold * (int(f_read[8]) - int(f_read[7]))):
						if sign == '-':
							operon.reverse()
						name = str(operon)
						if name in count_operon:
							count_operon[name] += 1
						else:
							count_operon[name] = 1
							operon_name[name] = operon
					elif len(operon) > 2:
						operon = operon[1:] #remove the first gene that is not the same sign as checked above
						if sign == '-':
							operon.reverse()
						name = str(operon)
						if name in count_operon:
							count_operon[name] += 1
						else:
							count_operon[name] = 1
							operon_name[name] = operon
					else:
						name = operon[1] #counting individual genes
						if name in count_operon:
							count_operon[name] += 1
						else:
							count_operon[name] = 1
				else:

					name = operon[0]
					if name in count_operon:
						count_operon[name] += 1
					else:
						count_operon[name] = 1
				prev = read[3]
				sign = read[5]
				operon = [read[9]]
				f_read = read
		for op in count_operon:
			if op in operon_name: #not a gene but an operon
				name_o = operon_name[op][0]+'-'+operon_name[op][-1]
				all_genes = str(operon_name[op])[1:-1]
				all_genes = all_genes.replace("'","")
				new.write('%s\t%s\t'%(name_o,all_genes))
				current_relative = list()
				for gene in operon_name[op]:
					if gene not in count_operon or count_operon[gene]==0:
						current_relative.append(1) #abundance is 1 because all the trasncription is in the form of operon not an individaul gene
					else:
						current_relative.append(count_operon[op]/count_operon[gene])

				new.write('%i\t'%(count_operon[op]))
				all_values = ''
				all_values = str(current_relative)[1:-1]
				new.write('%.2f\t' % (max(current_relative)))
				for value in current_relative[0:-1]:
					new.write('%.2f,'%(value))
				new.write('%.2f'%(current_relative[-1]))
				new.write('\n')
def reconstruct(geneBedFile,rawTranscriptFile,trIntesectFile,rawBedTrFile,rawBedTrFilesorted,outputFile):
	global config
	global log
	count = 0
	with open(rawTranscriptFile,'r') as old, open(rawBedTrFile,'w') as new:
		for line in old:
			count+=1
			read = line.strip().split()
			chr = read[0]
			start = int(read[1])
			end = int(read[2])
			both = [start,end]
			both.sort()
			nstart = both[0]
			nend = both[1]
			name = "Trasncript" +str(count) + "C"+ read[3]
			sign = read[4]
			new.write('%s\t%i\t%i\t%s\t60\t%s\n'%(chr,nstart,nend,name,sign))
	command="sort -k 1,1 -k2,2n "+rawBedTrFile +"> " +rawBedTrFilesorted
	process = subprocess.Popen(command, shell=True)
	process.wait()
	command=config["bedtools_path"]+" intersect -wao -a " +rawBedTrFilesorted +" -b " + geneBedFile + " > "+ trIntesectFile
	log.write("Running bedtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	with open(trIntesectFile,'r') as old, open(outputFile,'w') as new:
		line = old.readline()
		fread = line.strip().split('\t')
		prev = fread[3]
		sign = fread[5]
		start = int(fread[1])
		end = int(fread[2])
		length=end-start
		current = list()
		asCurrent = list()
		chromosome = fread[0]
		count = int(fread[3].split('C')[1])
		if fread[11] == sign:
			if (length > 0.2 *int(fread[8])-int(fread[7])):
				current.append(fread[9])
		else:
			if (length > 0.2 *int(fread[8])-int(fread[7])):
				asCurrent.append(fread[9])
		for line in old:
			read = line.strip().split()
			if read[3] == prev:					
				if read[11] == sign:
					if (length > 0.2 *int(read[8])-int(read[7])):
						current.append(read[9])
				else:
					if (length > 0.2 *int(read[8])-int(read[7])):
						asCurrent.append(read[9])
			else:
				if len(current) == 0:
					if len(asCurrent) ==0 or (len(asCurrent) == 1 and asCurrent[0] == '.'):
						new.write('%s\t%i\t%i\t%s\t%s\t%i\tNovel_Intergenic\n'%(chromosome,start,end,prev,sign,count))
					else:
						new.write('%s\t%i\t%i\t%s\t%s\t%i\tNovel_Antisesne\t%s\n'%(chromosome,start,end,prev,sign,count,asCurrent))
				else:
					new.write('%s\t%i\t%i\t%s\t%s\t%i\tOverlapping\t%s\n'%(chromosome,start,end,prev,sign,count,current))
				prev = read[3]
				sign = read[5]				
				current = list()
				asCurrent = list()
				chromosome = read[0]
				start = int(read[1])
				end = int(read[2])
				length=end-start
				count = int(read[3].split('C')[1])
				if read[11] == sign:
					if (length > 0.2 *int(read[8])-int(read[7])):
						current.append(read[9])
				else:
					if (length > 0.2 *int(read[8])-int(read[7])):
						asCurrent.append(read[9])
				


			






def seperateBedandGenomcov(bedfile,pos_bed,neg_bed,pos_cov,neg_cov,chrFile):
	global config
	global log
	# seperate pos and negative reads to calculate coverage seperately
	with open(bedfile, 'r') as old, open(pos_bed, 'w') as pos, open(neg_bed, 'w') as neg:
		for line in old:
			if line[-2] == '+':
				pos.write(line)
			else:
				neg.write(line)
	
	commandPos=config["bedtools_path"]+" "+config["genomecov_options"]+" -g " + chrFile+ " -i "+pos_bed+" > "+pos_cov
	log.write("Running bedtools with the following command: "+commandPos)
	process = subprocess.Popen(commandPos, shell=True)
	process.wait()
	commandNeg=config["bedtools_path"]+" "+config["genomecov_options"]+" -g " + chrFile+ " -i "+neg_bed+" > "+neg_cov
	log.write("Running bedtools with the following command: "+commandNeg)
	process = subprocess.Popen(commandNeg, shell=True)
	process.wait()
	log.write("Finished running bedtools")
def minimap(sourceFile, genomeFile, outputFile):
	global config
	global log
	command=config["minimap_path"]+" "+config["minimap_options"]+" "+genomeFile+" "+sourceFile+" > "+outputFile
	log.write("Running minimap with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running minimap")

def samToBam(sourceFile, outputFile, ChrFile):
	global config
	global log
	command=config["samtools_path"]+" "+config["samtobam_options"]+" "+sourceFile+" > "+outputFile
	log.write("Running samtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running samtools")
	#generating genome file (chromosomes and their lengths for later use)
	with open(sourceFile,'r') as old, open(ChrFile,'w') as new:
		for line in old:
			if line[0] == '@':
				read = line.strip().split()
				if read[0] == "@SQ":
					new.write('%s\t%i\n'%(read[1][3:],int(read[2][3:])))

def sortBam(sourceFile, outputFile):
	global config
	global log
	command=config["samtools_path"]+" "+config["sortbam_options"]+" "+sourceFile+" -o "+outputFile+" > "+outputFile
	log.write("Running samtools with the following command: "+command)
	process = subprocess.Popen(command, shell=True)
	process.wait()
	log.write("Finished running samtools")



def clustering(a, threshold):
	final = list() # clustering can produce more than TTS
	prev = a[0]
	curr = [a[0]]
	for i in a[1:]:
		if i - prev < threshold: #as long as the distance between reads are less than ~ bp, take them in the same cluster
			prev = i
			curr.append(i)
		else:
			if len(curr) > 0.2*len(a): #if the current cluster count is significant
				final.append((max(set(curr), key=curr.count), len(curr))) #then report the most frequent site as TSS
			prev = i
			curr = [i] #start a new cluster
	if len(curr) > 0.2 * len(a):
		final.append((max(set(curr), key=curr.count),len(curr)))
	return(final)




if __name__ == "__main__":
	main(sys.argv[1:])
