SYNOPSIS

This script reads a bam file and computes the reference base composition and substitution frequencies around 5p and 3p ends, including for di-nucleotides. Substitutions frequencies include all possible types; in addition CG->TG substitutions are considered separately. Reference base composition around alignment starts and ends can only be computed if the reference genome is provided. Input are sorted alignment files in BAM format produced with BWA.

USAGE
./substitution_patterns [-options] in.bam in2.bam ...

OPTIONS
-reference          provide reference genome in order to be able to compute reference base composition
                    outside of sequence alignments [optional; choose whole genome fasta file]
-minread            lower sequence length cutoff [default 30]
-maxread            upper sequence length cutoff [default 1000]
-quality            map quality cutoff [default 0]
-unpaired           do not exclude unpaired reads [off]
-graphic            produce graphical output
-xrange             define x-axis range for plotting [default 50]
-help               generates this help message
-report             define interval at which files are written [default 100_000]

OUTPUT
in.5p_refbase_composition.txt            reference base composition around 5p end (position 0 = terminal 5p base)
in.3p_refbase_composition.txt            reference base composition around 3p end (position 0 = terminal 3p base)
in.5p_substitutions.txt                  substitution frequencies at 5p end 
in.3p_substitutions.txt                  substitution frequencies at 3p end
in.Xp_dinucleotide_substitutions.txt     di-nucleotide substitution frequencies
in.Xp_dinucleotide_refbase_compositions  di-nucleotide reference base composition
in.plots.pdf                             graphical output
STDOUT                                   number of sequences going into analysis, number of sequences filtered out

DEPENDENCIES
This perl script has been successfully used with perl 5 (version 22). It requires 2 perl modules to be installed (File::Basename, Getopt::Long) as well as samtools (successfully used with version 1.3.1). If graphical output is turned on, gnuplot (version 5.0) and GPL Ghostscript (version 9.26) have to be installed.
