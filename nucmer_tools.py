import os, glob, sys
sys.path.append('/Users/like/scripts/tools/')
import fasta_tools, blast_tools

debug = False


#for a given id --> sequence dictionary, sort the ids according to size, returns an ordered list of contig-ids, plus a list of contigs that are smaller than minKb
def size_sort_contigs(cid2seq, minKb = 20):
	
	# first put contigs into a size --> contig dictionary
	size2contigs = {}
	for cid in id2seq.keys():	
		size = len(cid2seq[cid])
		if size2contigs.has_key(size): size2contigs[size].append(cid)
		else: size2contigs[size] = [cid]
		
	#then sort the contigsizes
	sizes = size2contigs.keys()
	sizes.sort()
	
	#sort() sorts from small to large, so we need to reverse this list
	sizes = sizes[::-1]
	
	sizeSorted_contigIDs = []
	small_contigIDs      = []
	for s in sizes:
		contigs = size2contigs[s]
		for c in contigs:
			sizeSorted_contigIDs.append(c)
			if s < minKb*1000: small_contigIDs.append(c)
				
	return sizeSorted_contigIDs, small_contigIDs



#When doing all x all genome alignments, we compare X versus Y but do not run MUMmer on Y vs X, to save time.
#Here we construct Y vs X coordinates-files from X vs Y data, because that is just a lot easier in downstream analyses (we don't have to constantly check what is the query and what is the reference)	
def complete_set_of_coords_files(list_of_genomes, dirname, extension):
	
	coords_files = set(glob.glob(dirname+'*'+extension))
	files_added  = []
	for x in range(len(list_of_genomes)):
		genomeX = list_of_genomes[x].split('/')[-1].split('.fasta')[0]
		for y in range(len(list_of_genomes)):
			genomeY = list_of_genomes[y].split('/')[-1].split('.fasta')[0]
			
			print genomeX, genomeY
			coords_fname_XY = dirname+genomeX+'.vs.'+genomeY+extension 
			if coords_fname_XY not in coords_files: # we don't have output for this comparison, we should have output for the reverse comparison
				coords_fname_YX = dirname+genomeY+'.vs.'+genomeX+extension
				print coords_fname_XY, 'does not exist, will infer contents based on', coords_fname_YX
				if not coords_fname_YX in coords_files: 
					print coords_fname_YX
					print "***ERROR***: can't find output for ", genomeY, 'vs', genomeX
					return 256
					
					
				lines = open(coords_fname_YX).readlines()
				Yfasta_fname, Xfasta_name = lines[0].strip().split() 
				#new_lines = lines[:5]
				
				#put new lines into a dictionary: Qscaffold --> start_hit_region in Qscaffold -> new line
				#so the output file can be ordered 
				#(coords files are ordered according to the reference fasta, here we swap query and reference, 
				#hence we also need to change the order) 
				scaffold2starts2new_lines = {}
				
				#file structure (header)
				#Qfasta_fname	Rfasta_fname
				#NUCMER
				#
				#[S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
				#=====================================================================================
				#  0	   1	2	   3	     4   5 	  6 	   7 	  8    9     10   11 12
				for line in lines[5:]: #skip header:
					data = line.strip().split()
					#print data
					s,e = map(int, data[3:5])
					
					new_line = ''
					if s<e:	new_line = data[3]+'\t'+data[4]+'\t|\t'+data[0]+'\t'+data[1]
					else:	new_line = data[4]+'\t'+data[3]+'\t|\t'+data[1]+'\t'+data[0]
					new_line += '\t|\t'+data[7]+'\t'+data[6]+'\t|\t'+data[9]+'\t|\t'+data[12]+'\t'+data[11]+'\n'
					
					Xscaffold = data[12]
					start = min([s,e])
					
					if scaffold2starts2new_lines.has_key(Xscaffold):
						if scaffold2starts2new_lines[Xscaffold].has_key(start):
							scaffold2starts2new_lines[Xscaffold][start].append(new_line)
						else:	scaffold2starts2new_lines[Xscaffold][start]=[new_line]
					else:
						scaffold2starts2new_lines[Xscaffold] = {}
						scaffold2starts2new_lines[Xscaffold][start]=[new_line]
						
				
				# get order of scaffolds from fasta file of genome X
				Xscaffold_list = []
				fasta_lines = open(Xfasta_name).readlines()
				for fline in fasta_lines:
					if fline[0] == '>': 
						#print fline.split()[0][1:]
						Xscaffold_list.append(fline[1:].split()[0].strip())
				
				print Xscaffold_list
				print
				print
				print scaffold2starts2new_lines.keys()
				# now create the new coords file
				outfile = open(coords_fname_XY, 'w')
				#write header
				outfile.write(Xfasta_name+'\t'+Yfasta_fname+'\n')
				for line in lines[1:5]:	outfile.write(line)
				for Xscaffold in Xscaffold_list: #use order of scaffolds from fasta
					if scaffold2starts2new_lines.has_key(Xscaffold):
						starts = scaffold2starts2new_lines[Xscaffold].keys()
						starts.sort()	#order start sites form small to large (so start at beginning of scaffold)
						for s in starts:
							for new_line in scaffold2starts2new_lines[Xscaffold][s]:
								outfile.write(new_line)
				outfile.close()
				
				files_added.append(coords_fname_XY)
				
	return files_added
							
			
#When doing all x all genome alignments, we compare X versus Y but do not run MUMmer on Y vs X, to save time.
#Here we construct Y vs X coordinates-files from X vs Y data, because that is just a lot easier in downstream analyses (we don't have to constantly check what is the query and what is the reference)	
def complete_set_of_coords_files_in_dir(dirname, extension):
	
	coords_files = set(glob.glob(dirname+'*'+extension))
	files_added  = []
	for coords_fname_YX in coords_files: 
		genomeY, genomeX = coords_fname_YX.split(dirname)[-1].split(extension)[0].split('.vs.')
		print genomeY, genomeX
		coords_fname_XY = dirname+genomeX+'.vs.'+genomeY+extension 
		if coords_fname_XY not in coords_files: # we don't have output for this comparison, we should have output for the reverse comparison
			print coords_fname_XY, 'does not exist, will infer contents based on', coords_fname_YX
			
			lines = open(coords_fname_YX).readlines()
			Yfasta_fname, Xfasta_name = lines[0].strip().split() 
			#new_lines = lines[:5]
			
			#put new lines into a dictionary: Qscaffold --> start_hit_region in Qscaffold -> new line
			#so the output file can be ordered 
			#(coords files are ordered according to the reference fasta, here we swap query and reference, 
			#hence we also need to change the order) 
			scaffold2starts2new_lines = {}
			
			#file structure (header)
			#Qfasta_fname	Rfasta_fname
			#NUCMER
			#
			#[S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
			#=====================================================================================
			#  0	   1	2	   3	     4   5 	  6 	   7 	  8    9     10   11 12
			for line in lines[5:]: #skip header:
				data = line.strip().split()
				#print data
				s,e = map(int, data[3:5])
				
				new_line = ''
				if s<e:	new_line = data[3]+'\t'+data[4]+'\t|\t'+data[0]+'\t'+data[1]
				else:	new_line = data[4]+'\t'+data[3]+'\t|\t'+data[1]+'\t'+data[0]
				new_line += '\t|\t'+data[7]+'\t'+data[6]+'\t|\t'+data[9]+'\t|\t'+data[12]+'\t'+data[11]+'\n'
				
				Xscaffold = data[12]
				start = min([s,e])
				
				if scaffold2starts2new_lines.has_key(Xscaffold):
					if scaffold2starts2new_lines[Xscaffold].has_key(start):
						scaffold2starts2new_lines[Xscaffold][start].append(new_line)
					else:	scaffold2starts2new_lines[Xscaffold][start]=[new_line]
				else:
					scaffold2starts2new_lines[Xscaffold] = {}
					scaffold2starts2new_lines[Xscaffold][start]=[new_line]
					
			
			# get order of scaffolds from fasta file of genome X
			Xscaffold_list = []
			fasta_lines = open(Xfasta_name).readlines()
			for fline in fasta_lines:
				if fline[0] == '>': 
					#print fline.split()[0][1:]
					Xscaffold_list.append(fline[1:].split()[0].strip())
			
			#print Xscaffold_list
			#print
			#print
			#print scaffold2starts2new_lines.keys()
			# now create the new coords file
			outfile = open(coords_fname_XY, 'w')
			#write header
			outfile.write(Xfasta_name+'\t'+Yfasta_fname+'\n')
			for line in lines[1:5]:	outfile.write(line)
			for Xscaffold in Xscaffold_list: #use order of scaffolds from fasta
				if scaffold2starts2new_lines.has_key(Xscaffold):
					starts = scaffold2starts2new_lines[Xscaffold].keys()
					starts.sort()	#order start sites form small to large (so start at beginning of scaffold)
					for s in starts:
						for new_line in scaffold2starts2new_lines[Xscaffold][s]:
							outfile.write(new_line)
			outfile.close()
			
			files_added.append(coords_fname_XY)
			
	return files_added			
			
##############################################
#
#Below are a number of functions that will determine in which order and direction 
#we have to place scaffolds in genome Q for it to best match genome R:
#i.e. how to get closest to a diagonal when plotting the whole genome alignment
#
##############################################

#for each scaffold/chr in the query, determine to which scaffold/chr in the reference it is best aligned, 
#based on the file with coordinates of hits between genome Q and genome R
#and the list of scaffolds in R
# for each scaffold/chr in the query, determine to which scaffold/chr in the reference it is best aligned.
def assign_scaffold2scaffold(coordsfilename, Rscaffolds, promer = False):
	
	# determine which scaffold belongs to which chromosome and in which direction
	scaffold2scaffold = {}
	lines = open(coordsfilename).readlines()[5:]
	Rscaffold2coordslines = {}
	for line in lines:
		data = line.split()
		Qscaffold = data[-1].strip()
		Rscaffold = data[-2].strip()
		#print data
		if not Rscaffold2coordslines.has_key(Rscaffold):
			Rscaffold2coordslines[Rscaffold] = [line]
		else:	Rscaffold2coordslines[Rscaffold].append(line)
		
		
		if not scaffold2scaffold.has_key(Qscaffold):
			scaffold2scaffold[Qscaffold] = {}
			for Rscaffold_tmp in Rscaffolds: scaffold2scaffold[Qscaffold][Rscaffold_tmp] = [[],[]] # two values, one for the + and one for the - strand
			
		#collect regions that correspond to alignments between the two scaffolds
		#print Qscaffold, 'assigned to', Rscaffold
		if int(data[0]) < int(data[1]):
			if int(data[3]) < int(data[4]): 
				scaffold2scaffold[Qscaffold][Rscaffold][0].append((int(data[0]), int(data[1])))
			else:	scaffold2scaffold[Qscaffold][Rscaffold][1].append((int(data[0]), int(data[1])))
		else:
			if int(data[3]) < int(data[4]): 
				scaffold2scaffold[Qscaffold][Rscaffold][1].append((int(data[1]), int(data[0])))
			else:	scaffold2scaffold[Qscaffold][Rscaffold][0].append((int(data[1]), int(data[0])))

	#Rscaffold2Qscaffolds = {}
	Qscaffold2direction  = {}
	Qscaffold2bestmatchingRscaffold = {}
	#for Rscaffold in Rscaffolds: Rscaffold2Qscaffolds[Rscaffold] = set([])
	
	for Qscaffold in scaffold2scaffold.keys():
		maxlength  = 0
		bestmatch  = ''
		maxlength0 = 0
		maxlength1 = 0
		#now check to which of the Reference scaffolds this Query scaffold has the longest alignment
		for Rscaffold in Rscaffolds:
			regions0 = blast_tools.merge_overlapping_regions(scaffold2scaffold[Qscaffold][Rscaffold][0], min_overlap = 1, minusstrand = False)
			regions1 = blast_tools.merge_overlapping_regions(scaffold2scaffold[Qscaffold][Rscaffold][1], min_overlap = 1, minusstrand = False)
			
			length0 = 0
			for (s,e) in regions0:
				length0 += e-s
			length1 = 0
			for (s,e) in regions1:
				length1 += e-s
				
			alignmentlength = length0 + length1
			
			#print Qscaffold, Rscaffold, alignmentlength
			if alignmentlength > maxlength:
				maxlength  = alignmentlength
				bestmatch  = Rscaffold
				maxlength0 = length0
				maxlength1 = length1
				
		if len(bestmatch)>0:
			
			# now determine the main direction of the alignments:
			strand = '+'
			if maxlength1 > maxlength0:
				strand = '-'
				
			#print Qscaffold, bestmatch, maxlength, strand
			Qscaffold2direction[Qscaffold] = strand
			Qscaffold2bestmatchingRscaffold[Qscaffold] = bestmatch
			#Rscaffold2Qscaffolds[bestmatch].add(Qscaffold)
	
	return Rscaffold2coordslines, scaffold2scaffold, Qscaffold2direction, Qscaffold2bestmatchingRscaffold
	


def get_order_and_direction_of_Qscaffolds__to_get_a_diagonal(coordsfile, Rscaffolds, outfilename=None, promer = False):
	
	# determine which scaffold belongs to which chromosome and in which direction
	Rscaffold2coordslines, scaffold2scaffold, Qscaffold2direction, Qscaffold2bestmatchingRscaffold = assign_scaffold2scaffold(coordsfile, Rscaffolds, promer = promer)
		
	mappedscaffolds = set([])
	
	out = ''
	Qscaffolds = []
	for Rscaffold in Rscaffolds:
		if debug: print Rscaffold
		
		if Rscaffold2coordslines.has_key(Rscaffold):
			lines = Rscaffold2coordslines[Rscaffold]
			for line in lines:
				data = line.split()
				Qscaffold = data[-1]
				if Qscaffold not in mappedscaffolds:
					if debug: print '\t', Qscaffold
					if Qscaffold2bestmatchingRscaffold.has_key(Qscaffold):
						if Rscaffold == Qscaffold2bestmatchingRscaffold[Qscaffold]:
							if debug: print 'best match!'
							out += Qscaffold+'\t'+Qscaffold2direction[Qscaffold]+'\n'
							mappedscaffolds.add(Qscaffold)
							Qscaffolds.append(Qscaffold)
							if debug: print Qscaffolds
			
						
	if outfilename != None:
		outfile = open(outfilename, 'w')
		outfile.write(out)
		outfile.close()
		
	return mappedscaffolds, Qscaffolds, Qscaffold2direction
	


if __name__ == "__main__":

	list_of_genomes = sys.argv[1].split(',')
	dirname = sys.argv[2]
	extension =  '.nucmer_maxmatch.coords'
	complete_set_of_coords_files(list_of_genomes, dirname, extension)




