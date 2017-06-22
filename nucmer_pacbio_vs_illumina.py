import glob, os, sys

Rname = ''
if len(sys.argv) > 1:
	Rname = sys.argv[1]

for pb_fasta in glob.glob('/Users/like/peter/WGAs/pacbio/'+Rname+'*.fasta'):
	pb_spid = pb_fasta.split('/')[-1].split('.fasta')[0]
	print pb_spid
	ifastas = glob.glob('/Users/like/peter/WGAs/fastas_genomes/'+pb_spid+'*.fasta')
	if len(ifastas) == 1:
		illumina_fasta = ifastas[0]
		cmnd = 'nice python /Users/like/scripts/genome_comparisons/run_nucmer__and__get_coords.py '+illumina_fasta+' /Users/like/peter/WGAs/compare_pacbio2illumina/ '+pb_fasta
		os.system(cmnd)

	cmnd = 'nice python /Users/like/scripts/genome_comparisons/run_nucmer__and__get_coords.py /Users/like/peter/WGAs/illumina/ /Users/like/peter/WGAs/nucmer_pacbio/ '+pb_fasta
	os.system(cmnd)

	