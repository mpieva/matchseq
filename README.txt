Author: Matthias Meyer
Date: May-15-2020

SYNOPSIS

This repository provides the following 3 perl scripts for reconstructing the structure of double-stranded DNA fragments from deeply sequenced single-stranded DNA libraries:

complementary_strands2.pl
alignment_distance.pl
substitution_patterns2.pl

Each script includes a help menu with further instructions and available options (accessible by executing the scripts without specifying an input file). 


SOFTWARE DEPENDENCIES

All script have been successfully used with perl 5 (version 22). They require 3 perl modules to be installed (File::Basename, Getopt::Long and List::Util), as well as samtools (successfully used with version 1.3.1). In addition, substitution_patterns2.pl requires gnuplot (version 5.0) and GPL Ghostscript (version 9.26) to be installed for graphical output (optional).


EXAMPLE 

In what follows we provide an example for how these scripts were used to determine the structure of DNA fragment ends, the length of gaps in the interior of DNA fragments and the base composition and substitution frequencies around fragment ends and gaps in the Vindija 33.19 Neandertal. The scripts and BAM files have to be located in the working directory or their path specified. 

We start with a BAM file (vindija33.19.bam, provided as part of the documentation) with sorted sequence alignments to the human reference genome (hg19). Sequences were pre-processed by adapter trimming and overlap-merging of paired reads (using leeHom, https://grenaud.github.io/leeHom/), removal of duplicates (using bam-rmdup, https://bitbucket.org/ustenzel/biohazard-tools) and removal of residual paired reads and sequences < 35 bp or with a map quality score < 25 (using ‘samtools view’).

/1/ Identification of overlapping sequence alignments on complementary strands (putatively representing double-stranded DNA fragments)

The complementary_strands2.pl script is executed as follows:

./complementary_strands2.pl -output 2 vindija33.19.bam

This script generates several output files, the most relevant of which are:

vindija33.19.comp_strands.hist.txt => tab-delimited table containing histograms of the distances between all 5' ends and the closest 3' end on the complementary strand (and vice versa).
vindija33.19.5p-blunt.bam => BAM file containing sequences where the 5’ end was inferred to be located in a blunt-end
vindija33.19.5p-overhang-1bp.bam => BAM file containing sequences where the 5’ end was inferred to be located in a 5’ overhang

Additional BAM files are generated for longer 5’ overhangs, 5’ recessed ends and sequences representing overhanging, blunt and recessed 3’ ends.

/2/ Identification of strands mapping in close proximity to the reference genome, in same strand orientation (putatively representing gaps in double-stranded DNA fragments)

Execute the alignment_distance.pl script as follows:

./alignment_distance.pl -output 2 vindija33.19.bam >vindija33.19.gaps.txt

This generates the following files:

vindija33.19.gaps.txt => tab delimited file containing a histogram of the shortest distances between the 3’ and 5’ ends of DNA strands mapped to the reference genome in the same orientation
vindija33.19.1st_dist2.bam => BAM file containing sequences with 3’ ends located 2 bp upstream of the 5’ end of another sequence mapped to the same strand (i.e. sequences preceding a 1 bp gap)
vindija33.19.1st_dist2.bam => BAM file containing sequences with 5’ ends located located 2 bp downstream of the 3’ end of another sequence mapped to the same strand (i.e. sequences preceded by a 1 bp gap)

/3/ Evaluating the position-dependent base composition and substitution frequencies in DNA strands assigned to various structural contexts.

Execute the substitution_frequencies2.pl script as follows, using the BAM files produced with the previous scripts as input. The reference genome (hg19) has to be provided in fasta format.

./substitution_patterns2.pl -graphical -reference hg19_whole_genome.fa *.bam

For each BAM file, a PDF file is created containing plots of the base composition and frequencies of all possible substitutions around alignment start and end points. The data underlying these plots are also provided in tab-delimited text format.