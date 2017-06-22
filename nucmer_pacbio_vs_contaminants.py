import glob, os, sys

Rname = ''
if len(sys.argv) > 1:
	Rname = sys.argv[1]

for pb_fasta in glob.glob('/Users/like/peter/WGAs/pacbio/'+Rname+'*.fasta'):
	pb_spid = pb_fasta.split('/')[-1].split('.fasta')[0]
	print pb_spid

	#'Usage: python run_nucmer__and__get_coords.py <dirname or fastafile> <outdirname> <ref_genome.fasta1,ref_genome.fasta2,..,ref_genome.fastan>'
	cmnd = 'nice python /Users/like/scripts/genome_comparisons/run_nucmer__and__get_coords.py /Users/like/peter/WGAs/contaminants/ /Users/like/peter/WGAs/contaminants/ '+pb_fasta
	os.system(cmnd)
	
	