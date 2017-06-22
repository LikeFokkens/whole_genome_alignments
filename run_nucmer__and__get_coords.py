PATH_MUMMER       = '/Applications/MUMmer3.23/'
PATH_PYTHON_TOOLS = '/Users/like/scripts/'
'''nice /Applications/MUMmer3.23/nucmer --maxmatch -p /Users/like/Shermineh/align_putative_centromeres/1361_lane7-index16.mapped2.Fol007.ordered4fol2.CenH3_regions__rseg.NUCMER /Users/like/Shermineh/align_putative_centromeres/1361_lane7-index16.mapped2.Fol007.ordered4fol2.CenH3_regions__rseg.fasta  /Users/like/Shermineh/align_putative_centromeres/1361_lane7-index16.mapped2.Fol007.ordered4fol2.CenH3_regions__rseg.fasta >& /Users/like/Shermineh/align_putative_centromeres/1361_lane7-index16.mapped2.Fol007.ordered4fol2.CenH3_regions__rseg.NUCMER.nucmer_log
'''

import glob, os, sys, time
sys.path.append(PATH_PYTHON_TOOLS)
import nucmer_tools

if len(sys.argv) == 1:
	print 'Usage: python run_nucmer__and__get_coords.py <dirname> <outdirname> <ref_genome.fasta1,ref_genome.fasta2,..,ref_genome.fastan>'
	print 'This script will compare all genomes in <dirname> to the <ref_genome.fasta>s'
	print 'If the <ref_genome> is ommitted, it will compare all against all'
	print 'All output will be saved in <outdirname>/deltafiles, <outdirname>/coords, log files in <outdirname>/logs'
	print 'fastas should end with ".fasta"'

	sys.exit()
	
dirname    = sys.argv[1]
fastafiles = [sys.argv[1]]							#
if os.path.isdir(dirname):							#  added 2015 11 24
	fastafiles = glob.glob(dirname + '/*.fasta')	#
outdirname = sys.argv[2]

if not os.path.exists(outdirname+'/deltafiles'): os.mkdir(outdirname+'/deltafiles')
if not os.path.exists(outdirname+'/coords'): os.mkdir(outdirname+'/coords')
if not os.path.exists(outdirname+'/logs'): os.mkdir(outdirname+'/logs')

for Qi, Qfasta in enumerate(fastafiles):
	Qname = Qfasta.split('/')[-1].split('.fa')[0]
	Rfastas = []
	
	if len(sys.argv)>3:
		if os.path.isdir(sys.argv[3]):							#  added 2015 11 24
			Rfastas = glob.glob(sys.argv[3] + '/*.fasta')
		else:
			Rfastas = sys.argv[3].split(',')
	else:
		Rfastas = fastafiles[Qi:]
	
	for Rfasta in Rfastas:
		Rname = Rfasta.split('/')[-1].split('.fa')[0]
		nucmer_out_prefix    = outdirname+'/deltafiles/'+Rname+'.vs.'+Qname+'.nucmer_maxmatch'
		nucmer_out_prefixREV = outdirname+'/deltafiles/'+Qname+'.vs.'+Rname+'.nucmer_maxmatch'
		coords_out_prefix    = outdirname+'/coords/'+Rname+'.vs.'+Qname+'.nucmer_maxmatch'
		coords_out_prefixREV = outdirname+'/coords/'+Qname+'.vs.'+Rname+'.nucmer_maxmatch'

		nucmer_cmnd = 'nice '+PATH_MUMMER + 'nucmer --maxmatch -p '+nucmer_out_prefix+' '+Rfasta+' '+Qfasta+' >& '+outdirname+'/logs/'+Rname+'.vs.'+Qname+'.nucmer_maxmatch.log'
		coords_cmnd = 'nice '+PATH_MUMMER + 'show-coords -r '+nucmer_out_prefix+'.delta > '+coords_out_prefix+'.coords'
			
		if not os.path.exists(nucmer_out_prefix+'.delta') and not os.path.exists(nucmer_out_prefixREV+'.delta'):
			nucmerres = os.system(nucmer_cmnd)
		else: 
			print nucmer_out_prefix+'.delta or', nucmer_out_prefixREV+".delta exists, I will not execute [", nucmer_cmnd, "]"
			nucmerres = 'not executed'

		#time.sleep(3)
		if not os.path.exists(coords_out_prefix+'.coords') and not os.path.exists(coords_out_prefixREV+'.coords'):
			coordsres = os.system(coords_cmnd)
		else: 
			print coords_out_prefix+'.coords', "exists, I will not execute [", coords_cmnd, "]"
			coordsres = 'not executed'
		#print nucmer_cmnd
		print Rname, Qname, 'nucmer:', nucmerres,', coords:', coordsres
	


c_extension = '.nucmer_maxmatch.coords'
c_dirname   = outdirname+'/coords/'

if len(sys.argv) == 3:	
	inferred_files = nucmer_tools.complete_set_of_coords_files(fastafiles, c_dirname, c_extension)
else:
	inferred_files = nucmer_tools.complete_set_of_coords_files_in_dir(c_dirname, c_extension)

print len(inferred_files), 'coords-files added to complete all X all'
#just checking:
for fname_xy in inferred_files:
	len_xy = len(open(fname_xy).readlines())
	
	x,y = fname_xy.split('/')[-1].split(c_extension)[0].split('.vs.')
	fname_yx = c_dirname+y+'.vs.'+x+c_extension
	len_yx = len(open(fname_yx).readlines())
	
	print fname_xy, len_xy, fname_yx, len_yx
		
