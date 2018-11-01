# whole_genome_alignments

Contains the pipeline I used to create whole genome alignments between Fusarium oxysporum genomes. This is in 2014, using version of MUMmer (MUMMer3.23) and python (2.6 or 2.7) that were current at that time.

This is placed here as a surplus for the Materials and Methods in the paper [insert ref once published :-)], not as working standalone code or anything. Maybe later.

I used presence_absence_plots in genome-wide_plots to visualize results.

This code depends on 'nucmer_tools.py' that is in the tools repository. 
If you put this 'tools' dir in $PYTHONPATH it should work, or else you can use the following structure

path_to_your_core/tools
path_to_your_core/whole_genome_alignments

this should also work.







