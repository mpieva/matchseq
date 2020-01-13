#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;

my $qual_cutoff = 0;
my $minread = 30;
my $direction = 1;
my $range = 100;
my ($output, $side);

GetOptions ("quality=s" => \$qual_cutoff,
	    "minread=i" => \$minread,
            "range=i" => \$range,
            "output=i" => \$output);

&help unless scalar @ARGV == 1;
&help unless $ARGV[0] =~ /\.bam$/;
if ($output) {
    my $outname = basename($ARGV[0]);
    $outname =~ s/\.bam//;
    open WRITE5, "| samtools view -S -b - >$outname.1st_dist$output.bam" or die "oculd not write file $! \n";
    open WRITE3, "| samtools view -S -b - >$outname.2nd_dist$output.bam" or die "could not write file $! \n";
}

open BAM, "samtools view -hX $ARGV[0] |" or die "could not open file $!\n";

my (%dist_hist, %gap_base, %gap_before, %gap_after);
my ($counter, $interval) = (0, 0);
my $prev_chr = "";
my $prev_end_plus = 0;
my $prev_end_minus = 0;
my ($prev_line_plus, $prev_line_minus);
while (my $line = <BAM>) {
    $counter++; $interval++;
    if ($line =~ /^\@/) {
	if ($output) {
	    print WRITE5 $line;
	    print WRITE3 $line;
	}
	next;
    }
    if ($interval == 100_000) {
	print STDERR "\r$counter";
	$interval = 0;
    }
    chomp $line;
    my @fields = split /\t/, $line;
    my ($header, $flag, $chr, $start, $mapqual, $cigar, $junk1, $junk2, $junk3, $sequence) = split /\t/, $line;
    my $length = length($sequence);
    next unless $chr =~ /^(\d{1,2}|X)$/;
    next if $length < $minread;
    next if $mapqual < $qual_cutoff;
    my $end = aln_end($start, $cigar, $length);
    if ($flag =~ /r/) {
	if ($chr eq $prev_chr && $prev_end_minus > 0 && $prev_end_minus < $start) {
	    my $distance = $start - $prev_end_minus;
	    $dist_hist{$distance}++;
	    if ($output && $distance == $output) {
	    }
	}
	$prev_end_minus = $end;
	$prev_line_minus = $line;
    } else {
	if ($chr eq $prev_chr && $prev_end_plus > 0 && $prev_end_plus < $start) {
	    my $distance = $start - $prev_end_plus;
	    $dist_hist{$distance}++;
	    if ($output && $distance == $output) {
		print WRITE3 "$line\n";
		print WRITE5 "$prev_line_plus\n" if $prev_line_plus;
	    }
	}
	$prev_end_plus = $end;
	$prev_line_plus = $line;
    }
    $prev_chr = $chr;
}
foreach my $dist (1 .. $range) {
    print "$dist\t", $dist_hist{$dist} || 0, "\n";
}

sub aln_end {
    my ($start, $cigar, $length) = @_;
    my $end;
    if ($cigar =~ /[ID]/) {
        my ($insertions, $deletions, $matches) = (0,0,0);
        my $cigar_backup = $cigar;
        while ($cigar =~ s/(\d+)([MID])//) {
            my ($num, $type) = ($1, $2);
            $insertions += $num if $type eq "I";
            $deletions += $num if $type eq "D";
            $matches += $num if $type eq "M";
        }
        $end = $start + $matches + $deletions - 1;
    } else {
        $end = $start + $length - 1;
    }
    return $end;
}

sub help {
print STDERR "

This script computes a histogram of the distance between the 3' end of each sequence and the closest 5' end of another sequence on the same strand. A 1-bp gap corresponds to a distance of 2. Input are sorted alignment files in BAM format produced by BWA.

              --------------->          --------------->           ----------->
              5              3          5              3           5          3
                        ---------------------->
                        5                     3

                             |<-------->|              |<--------->|
                        

[usage]
alignment_distance.pl

[options]
-quality      map quality cutoff [0]
-minread      minimal length cutoff [35]
-range        range of values to report
-output       write out sequences that have another sequence in distance X from their 5' or 3' end respectively. 

[output]
infile.1st_distX.bam   first (left) sequence whose 3' end is in distance X to the 5' end of the following sequence
infile.2nd_distX.bam   second (right) sequences whose 5' end is in distance X to the 3' of the sequence to the left
screen                 alignment distance histogram
";
exit;
}
