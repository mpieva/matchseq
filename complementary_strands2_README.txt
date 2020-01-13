SYNOPSIS

This script computes a histogram of the distance between the 5' and 3' end of each strand to 3' and 5' end of the closest complementary strand. It also computes the closest distance between 5'/3' ends that are on the same strand for comparison. If desired, strands that are in distance X to the closest 3' or 5' end of the complementary strand can be written to separate bam files. Only alignments to the autosomes or chromosome X are analyzed. Input are sorted alignment files in BAM format created with BWA. 

CAUTION: For mate pairs, only the first read of the pair is currently written to bam output. Histograms are computed correctly. 

USAGE
./complementary_strands2.pl [-options] in.bam

[options]
-quality      map quality cutoff [0]
-minread      minimal length cutoff [30]
-range        range of values to report [100]
-output       write out sequences that have another sequence in distance X from their 5' or 3' end respectively. 
-control      write out sequences that have another sequence in distance X from their 5' or 3' end respectively IN THE SAME ORIENTATION
-mt           analyze mtDNA (drop mapping restrictions to X and autosomes)

OUTPUT
-in.hist.txt  table with histogram of closest distances, in number of observations and fraction of sequences (>redirect to file)
-in.ends.txt  table showing the observed combinations between 5' and 3' molecule structures
-in.matching_hist.txt table providing the fraction of sequences with a maximum distance of X between complementary strands, binned by sequence length 
-bam files    in.{5p/3p}-{overhang/recessed/blunt}-[X]bp.bam

DEPENDENCIES
This perl script has been successfully used with perl 5 (version 22). It requires 3 perl modules to be installed (File::Basename, Getopt::Long and List::Util) as well as samtools (successfully used with version 1.3.1). 
