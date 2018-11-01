PATH_MUMMER       = '/Applications/MUMmer3.23/'

'''nice /Applications/MUMmer3.23/nucmer --maxmatch -p /Users/like/Shermineh/align_putative_centromeres/1361_lane7-index16.mapped2.Fol007.ordered4fol2.CenH3_regions__rseg.NUCMER 
/Users/like/Shermineh/align_putative_centromeres/1361_lane7-index16.mapped2.Fol007.ordered4fol2.CenH3_regions__rseg.fasta  /Users/like/Shermineh/align_putative_centromeres/1361_lane7-index16.mapped2.Fol007.ordered4fol2.CenH3_regions__rseg.fasta >& /Users/like/Shermineh/align_putative_centromeres/1361_lane7-index16.mapped2.Fol007.ordered4fol2.CenH3_regions__rseg.NUCMER.nucmer_log
'''
import glob, os, sys, time
# import tools sets:
# get path of this script, assume tools are in directory 'tools/' in the parent directory of this code
parent_dir  = os.path.dirname(os.path.abspath(os.path.dirname(sys.argv[0])))
PATH_PYTHON_TOOLS = parent_dir+'/tools/'


sys.path.append(PATH_PYTHON_TOOLS)
import nucmer_tools

if len(sys.argv) == 1:
	print('Usage: python run_promer__and__get_coords.py <dirname> <outdirname> <ref_genome.fasta1,ref_genome.fasta2,..,ref_genome.fastan>')
	print('This script will compare all genomes in <dirname> to the <ref_genome.fasta>s')
	print('If the <ref_genome> is ommitted, it will compare all against all')
	print('All output will be saved in <outdirname>/deltafiles, <outdirname>/coords, log files in <outdirname>/logs')
	print('fastas should end with ".fasta" or ".fa"')

	sys.exit()
	
dirname    = sys.argv[1]
outdirname = sys.argv[2]
fastafiles = glob.glob(dirname + '/*.fasta') + glob.glob(dirname + '/*.fa')


if not os.path.exists(outdirname+'/deltafiles'): os.mkdir(outdirname+'/deltafiles')
if not os.path.exists(outdirname+'/coords'): os.mkdir(outdirname+'/coords')
if not os.path.exists(outdirname+'/logs'): os.mkdir(outdirname+'/logs')

for Qi, Qfasta in enumerate(fastafiles):

	Qname = Qfasta.split('/')[-1].split('.fa')[0]
	Rfastas = []
	
	if len(sys.argv)>3:
		Rfastas = sys.argv[3].split(',')
	else:
		Rfastas = fastafiles[Qi:]
	

	for Rfasta in Rfastas:
		if not (os.path.exists(Rfasta)): 
			print ("Can not find ", Rfasta, "\nPlease check your input, did you add the correct path?")
			sys.exit()
			
		Rname = Rfasta.split('/')[-1].split('.fa')[0]
		promer_out_prefix    = outdirname+'/deltafiles/'+Rname+'.vs.'+Qname+'.promer_maxmatch'
		promer_out_prefixREV = outdirname+'/deltafiles/'+Qname+'.vs.'+Rname+'.promer_maxmatch'
		coords_out_prefix    = outdirname+'/coords/'+Rname+'.vs.'+Qname+'.promer_maxmatch'
		coords_out_prefixREV = outdirname+'/coords/'+Qname+'.vs.'+Rname+'.promer_maxmatch'

		promer_cmnd = 'nice '+PATH_MUMMER + 'promer --maxmatch -p '+promer_out_prefix+' '+Rfasta+' '+Qfasta+' >& '+outdirname+'/logs/'+Rname+'.vs.'+Qname+'.promer_maxmatch.log'
		coords_cmnd = 'nice '+PATH_MUMMER + 'show-coords -r '+promer_out_prefix+'.delta > '+coords_out_prefix+'.coords'	
		
		if not os.path.exists(promer_out_prefix+'.delta') and not os.path.exists(promer_out_prefixREV+'.delta'):
			promerres = os.system(promer_cmnd)
		else: 
			print(promer_out_prefix+'.delta or', promer_out_prefixREV+".delta exists, I will not execute [", promer_cmnd, "]")
			promerres = 'not executed'

		time.sleep(3)
		if not os.path.exists(coords_out_prefix+'.coords') and not os.path.exists(coords_out_prefixREV+'.coords'):
			coordsres = os.system(coords_cmnd)
		else: 
			print(coords_out_prefix+'.coords', "exists, I will not execute [", coords_cmnd, "]")
			coordsres = 'not executed'
		
		print(Rname, Qname, 'promer:', promerres,', coords:', coordsres)
	

if len(sys.argv) == 3:	
	c_extension = '.promer_maxmatch.coords'
	c_dirname   = outdirname+'/coords/'
	inferred_files = promer_tools.complete_set_of_coords_files(fastafiles, c_dirname, c_extension)
	
	print(len(inferred_files), 'coords-files added to complete all X all')
	#just checking:
	for fname_xy in inferred_files:
		len_xy = len(open(fname_xy).readlines())
		
		x,y = fname_xy.split('/')[-1].split(c_extension)[0].split('.vs.')
		fname_yx = c_dirname+y+'.vs.'+x+c_extension
		len_yx = len(open(fname_yx).readlines())
		
		print(fname_xy, len_xy, fname_yx, len_yx)
		

