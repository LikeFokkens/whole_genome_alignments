# whole_genome_alignments

Contains the pipeline I used to create whole genome alignments between Fusarium oxysporum genomes. This is in 2014, using version of MUMmer (MUMMer3.23) and python (2.6 or 2.7) that were current at that time.

This is placed here as a surplus for the Materials and Methods in the paper [insert ref once published :-)], not as working standalone code or anything. Maybe later.

I used presence_absence_plots in genome-wide_plots to visualize results.

This code depends on 'nucmer_tools.py' that is in the tools repository. 
If you put this 'tools' dir in $PYTHONPATH it should work, or else you can use the following structure. 

`path_to_your_core/tools`  
`path_to_your_core/whole_genome_alignments`  
  
this should also work.  
  
This code is a bit old, before I started using argument parser on a regular basis. If you want to know how to use it, just type the name of the file without arguments and you get some explanantion. E.g. `run_nucmer__and__get_coords.py` will print:  
` Usage: python run_nucmer__and__get_coords.py <dirname> <outdirname> <ref_genome.fasta1,ref_genome.fasta2,..,ref_genome.fastan>`    
`This script will compare all genomes in <dirname> to the <ref_genome.fasta>.`  
`If the <ref_genome> is ommitted, it will compare all against all.`  
`All output will be saved in <outdirname>/deltafiles, <outdirname>/coords, log files in <outdirname>/logs`  
`fastas should end with ".fasta" or ".fa"`. 





