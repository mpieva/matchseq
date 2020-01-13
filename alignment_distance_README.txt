SYNPOSIS

This script computes a histogram of the distance between the 3' end of each sequence and the closest 5' end of another sequence on the same strand. A 1-bp gap corresponds to a distance of 2. Input are sorted alignment files in BAM format produced by BWA.

USAGE
alignment_distance.pl

OPTIONS
-quality      map quality cutoff [0]
-minread      minimal length cutoff [35]
-range        range of values to report
-output       write out sequences that have another sequence in distance X from their 5' or 3' end respectively. 

OUTPUT
infile.1st_distX.bam   first (left) sequence whose 3' end is in distance X to the 5' end of the following sequence
infile.2nd_distX.bam   second (right) sequences whose 5' end is in distance X to the 3' of the sequence to the left
screen                 alignment distance histogram

DEPENDENCIES
This perl script has been successfully used with perl 5 (version 22). It requires 3 perl modules to be installed (File::Basename and Getopt::Long) as well as samtools (successfully used with version 1.3.1). 
