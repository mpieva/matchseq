#!/usr/bin/perl -w
use strict "vars";
use strict "subs";
use File::Basename;
use Getopt::Long;
use List::Util 'shuffle';

my $qual_cutoff = 0;
my $minread = 30;
my $direction = 1;
my $range = 100;
my $mt;
my $output;
my $control;


GetOptions ("quality=s" => \$qual_cutoff,
	    "minread=i" => \$minread,
            "range=i" => \$range,
	    "mt" => \$mt,
	    "control" => \$control,
            "output=i" => \$output);

&help unless scalar @ARGV == 1;
&help unless $ARGV[0] =~ /\.bam$/;

my (@seqs_to_analyze, @all_seqs, $bam_header, $counter, $interval, %dist_hist, %structure_hist, %match_hist, %filehandles);
my ($last_ms, $last_me, $last_ps, $last_pe, $last_chr) = (0, 0, 0, 0);
my ($first_ms, $first_me, $first_ps, $first_pe);

my $outname = basename($ARGV[0]);
$outname =~ s/\.bam//;
open BAM, "samtools view -h $ARGV[0] |" or die "could not open file\n";

&add_sequence;
while (scalar @seqs_to_analyze > 0) {
    my $current_seq = shift @seqs_to_analyze;
    my ($orientation, $chr, $start, $end, $line, $seqlength) = @{$current_seq};
    $counter++, $interval++;
    if ($interval == 10_000) {
        printf STDERR "\r$counter %20s", "processing chromosome $chr    ";
        $interval = 0;
    }
    while ("remove sequences from previous chromosome") {
	my $first_seq = shift @all_seqs;
	my ($seq_orientation, $seq_chr, $seq_start, $seq_end, $seq_line) = @{$first_seq};
	if ($seq_chr ne $chr) {
	    ($last_ps, $last_ms) = (undef, undef);
	    next;
	}
	unshift @all_seqs, $first_seq;
	last;
    }
    while ("keep adding sequences until enough") {
	if ($last_chr ne $chr) {
	    last;
	}
	if ($last_ps && $last_ms && $last_ps > $end && $last_ms > $end) {
	    last;
	}
	&add_sequence || last;
    }
    my @all_seqs_cleaned;
    my ($left_dist_ps, $left_dist_pe, $left_dist_ms, $left_dist_me);
    my ($first_pe, $first_me);
    while ("searching left and eliminating unused sequences") {
	my $next_seq = pop @all_seqs || last;
	unshift @all_seqs_cleaned, $next_seq;
	my ($seq_orientation, $seq_chr, $seq_start, $seq_end, $seq_line) = @{$next_seq};
	$first_pe = $seq_end if $seq_orientation eq "plus";
	$first_me = $seq_end if $seq_orientation eq "minus";
	next if $seq_line eq $line; #do not compare sequences to themselves
	next if $seq_chr ne $chr;
	if ($seq_start - $start <= 0) {
	    if ($seq_orientation eq "plus") {
		$left_dist_ps = $seq_start - $start unless defined $left_dist_ps;
	    } else {
		$left_dist_ms = $seq_start - $start unless defined $left_dist_ms;
 	    }
	} 
	if ($seq_end - $end <= 0) {
	    if ($seq_orientation eq "plus") {
		$left_dist_pe = $seq_end - $end unless defined $left_dist_pe;
	    } else {
		$left_dist_me = $seq_end - $end unless defined $left_dist_me;
	    }
	}
	if (defined $left_dist_ps && defined $left_dist_ms && defined $left_dist_pe && defined $left_dist_me) {
	    last if (defined $first_pe && defined $first_me && $first_pe < $start && $first_me < $start);
	}
    }
    @all_seqs = @all_seqs_cleaned;
    my ($right_dist_ps, $right_dist_pe, $right_dist_ms, $right_dist_me);
    while ("searching right sequences") {
IT:	foreach my $next_seq (@all_seqs) {
	    my ($seq_orientation, $seq_chr, $seq_start, $seq_end, $seq_line) = @{$next_seq};
	    next if $seq_line eq $line;
	    next if $seq_chr ne $chr;
	    if ($seq_start - $start > 0) {
		if ($seq_orientation eq "plus") {
		    $right_dist_ps = $seq_start - $start unless $right_dist_ps;
		} else {
		    $right_dist_ms = $seq_start - $start unless $right_dist_ms;
		}
	    }
	    if ($seq_end - $end > 0) {
		if ($seq_orientation eq "plus") {
		    $right_dist_pe = $seq_end - $end unless $right_dist_pe;
		} else {
		    $right_dist_me = $seq_end - $end unless $right_dist_me;
		}
	    }
	    if (defined $right_dist_ps && defined $right_dist_ms && defined $right_dist_pe && defined $right_dist_me) {
		last IT;
	    }
	}
	last;
    }
    
    my ($dist_5, $dist_3, $dist_5_control, $dist_3_control);
    if ($orientation eq "plus") {
	$dist_5 = &store_distance($left_dist_ms, $right_dist_ms, "dist_ps_ms", "dist_5");
	$dist_3 = &store_distance($left_dist_me, $right_dist_me, "dist_pe_me", "dist_3");
	$dist_5_control = &store_distance($left_dist_ps, $right_dist_ps, "dist_ps_ps", "dist_5_control");
	$dist_3_control = &store_distance($left_dist_pe, $right_dist_pe, "dist_pe_pe", "dist_3_control");
    } elsif ($orientation eq "minus") {
	$dist_3 = &store_distance($left_dist_ps, $right_dist_ps, "dist_ms_ps", "dist_3", "NEG");
	$dist_5 = &store_distance($left_dist_pe, $right_dist_pe, "dist_me_pe", "dist_5", "NEG");
	$dist_3_control = &store_distance($left_dist_ms, $right_dist_ms, "dist_ms_ms", "dist_3_control", "NEG");
	$dist_5_control = &store_distance($left_dist_me, $right_dist_me, "dist_me_me", "dist_5_control", "NEG");
    }
    if (defined $output && abs($dist_5) <= $output) {
	my $string = "5p-";
	$string .= "blunt" if $dist_5 == 0;
	$string .= "overhang-${dist_5}bp" if $dist_5 > 0;
	$string .= "recessed-" . abs($dist_5) . "bp" if $dist_5 < 0;
	my $fhdist = "fh_dist5_$dist_5";
	unless (exists $filehandles{$fhdist}) {
	    open ($fhdist, "| samtools view -bS - >$outname.$string.bam") or die "could not write file\n";
	    print $fhdist $bam_header;
	    print $fhdist $line;
	    $filehandles{$fhdist}++;
	} else {
	    print $fhdist $line;
	}
    }

    if (defined $output && abs($dist_3) <= $output) {
	my $string = "3p-";
	$string .= "blunt" if $dist_3 == 0;
	$string .= "overhang-" . abs($dist_3) . "bp" if $dist_3 < 0;
	$string .= "recessed-" . abs($dist_3) . "bp" if $dist_3 > 0;
	my $fhdist = "fh_dist3_$dist_3";
	unless (exists $filehandles{$fhdist}) {
	    open ($fhdist, "| samtools view -bS - >$outname.$string.bam") or die "could not write file\n";
	    print $fhdist $bam_header;
	    print $fhdist $line;
	    $filehandles{$fhdist}++;
	} else {
	    print $fhdist $line;
	}
    }

    if ($control && defined $output && abs($dist_5_control) <= $output) {
	my $string = "control_5p-";
	$string .= "blunt" if $dist_5_control == 0;
	$string .= "overhang-${dist_5_control}bp" if $dist_5_control > 0;
	$string .= "recessed-" . abs($dist_5_control) . "bp" if $dist_5_control < 0;
	my $fhdist = "fh_dist_5_control_$dist_5_control";
	unless (exists $filehandles{$fhdist}) {
	    open ($fhdist, "| samtools view -bS - >$outname.$string.bam") or die "could not write file\n";
	    print $fhdist $bam_header;
	    print $fhdist $line;
	    $filehandles{$fhdist}++;
	} else {
	    print $fhdist $line;
	}
    }
    
    if ($control && defined $output && abs($dist_3_control) <= $output) {
	my $string = "control_3p-";
	$string .= "blunt" if $dist_3_control == 0;
	$string .= "overhang-" . abs($dist_3_control) . "bp" if $dist_3_control < 0;
	$string .= "recessed-" . abs($dist_3_control) . "bp" if $dist_3_control > 0;
	my $fhdist = "fh_dist3_control_$dist_3_control";
	unless (exists $filehandles{$fhdist}) {
	    open ($fhdist, "| samtools view -bS - >$outname.$string.bam") or die "could not write file\n";
	    print $fhdist $bam_header;
	    print $fhdist $line;
	    $filehandles{$fhdist}++;
	} else {
	    print $fhdist $line;
	}
    }
    $structure_hist{"$dist_5,$dist_3"}++;
    my $min_dist = abs($dist_3);
    $min_dist = abs($dist_5) if abs($dist_5) < abs ($dist_3);
    foreach ($min_dist .. 100) {
	$match_hist{$seqlength}{$_}++;
    }
    $match_hist{$seqlength}{"obs"}++;
}

open HIST, ">$outname.comp_strands.hist.txt" or die "could not write histogram\n";
print HIST join ("\t", qw /Dist 5p_end %5p_end 3p_end %3p_end 5p_end_control %5p_end_control 3p_end 3p_end_control/), "\n";
foreach my $num (-$range..$range) {
    print HIST "$num";
    foreach my $type (qw/dist_5 dist_3 dist_5_control dist_3_control/) {
	print HIST "\t", $dist_hist{$type}{$num} || 0;
	my $percent;
	my $obs = $dist_hist{$type}{"obs"};
	if ($obs) {
	    my $count = $dist_hist{$type}{$num} || 0;
	    $percent = sprintf "%.3f", $count / $obs * 100;
	} else {
	    $percent = 0;
	}
	print HIST "\t$percent";
    }
    print HIST "\n";
}

open STRUC, ">$outname.ends.txt" or die "could not write ends\n";
print STRUC "|5p|\-3p-\t", join ("\t", -100 .. 100), "\n";
foreach my $dist5 (-100 .. 100) {
    print STRUC $dist5;
    foreach my $dist3 (-100 .. 100) {
	print STRUC "\t", $structure_hist{"$dist5,$dist3"} || 0;
    }
    print STRUC "\n";
}

open MATCH, ">$outname.matching_hist.txt" or die "could not write matching stats\n";
print MATCH join ("\t", "Length", 0 .. 30), "\n";
foreach my $seqlength (30 .. 100) {
    print MATCH $seqlength;
    foreach my $min_dist (0 .. 50) {
	my $num = $match_hist{$seqlength}{$min_dist} || 0;
	my $obs = $match_hist{$seqlength}{"obs"} || 0;
	if ($obs) {
	    my $fraction = sprintf ("%.5f", $num / $obs);
	    print MATCH "\t$fraction";
	} else {
	    print MATCH "\tNA";
	}
    }
    print MATCH "\n";
}

sub store_distance {
    my $distance;
    my ($left_dist, $right_dist, $key1, $key2, $negative) = @_;
    if (defined $left_dist && defined $right_dist && abs($left_dist) == abs($right_dist)) {
	my @tmp = shuffle ($left_dist, $right_dist);
	$distance = shift @tmp;
    } elsif (defined $left_dist && defined $right_dist) {
	$distance = $left_dist;
	$distance = $right_dist if abs($right_dist) < abs($left_dist);
    } elsif (defined $left_dist || defined $right_dist) {
	$distance = $left_dist;
	$distance = $right_dist if defined $right_dist;
    } else { 
	warn "OH: $left_dist, $right_dist\n";
    }
    $distance = -$distance if $negative;
    $dist_hist{$key1}{$distance}++;
    $dist_hist{$key1}{"obs"}++;
    $dist_hist{$key2}{$distance}++;
    $dist_hist{$key2}{"obs"}++;
    return $distance;
}

sub add_sequence {
    while (1) {
	my $line = <BAM> || return 0;
	if ($line =~ /^\@/) {
	    $bam_header .= "$line";
	    next;
	}
	my ($header, $flag, $chr, $start, $mapqual, $cigar, $junk1, $mate_start, $tlen, $sequence) = split /\t/, $line;
	$flag = decode_flag($flag);
	my $length = length($sequence);
	unless ($mt) {
	    next unless $chr =~ /^(chr)?(\d{1,2}|X)$/;
	}
	next if $length < $minread;
	next if $mapqual < $qual_cutoff;
	next if $flag =~ /u/i;
	my $end;
	if ($flag =~ /p/) {
	    next unless $flag =~ /1/; ##only ever look at the first read in a pair
	    next if $tlen > 2_000;
	    next if $tlen < 30;
	    unless ($flag =~ /r/) { ##read 1 on plus strand
#		die "unexpected paired read orientation\n:FLAG:$flag LINE:$line\n" if $tlen < 0;
		next if $tlen < 0;
		$end = $start + $tlen - 1;
	    } else { ##read 1 on minus strand
		next if $tlen > 0;
#		warn "unexpected reverse read orientation:\nFLAG:$flag LINE:$line" if $tlen > 0;
		$start = $mate_start;
		$end = $mate_start + abs($tlen) - 1;
	    }
	} else { #single reads
	    $end = &aln_end($start, $cigar, $length);
	}
	$last_chr = $chr;
	my $orientation;
	if ($flag =~ /r/) {
	    $orientation = "minus";
	    $last_ms = $start;
	} else {
	    $orientation = "plus";
	    $last_ps = $start;
	}
	my @seq = ($orientation, $chr, $start, $end, $line, $length);
	push (@seqs_to_analyze, \@seq);
	push (@all_seqs, \@seq);
	return 1;
    }
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
 
sub decode_flag {
    my $hex = shift(@_);
    my $bin = sprintf "%b", $hex;
    my $literal = "qdfs21RrUuPp";
    my $string;
    foreach my $pos (-length($bin) .. -1) {
	$string .= substr $literal, $pos, 1 if (substr $bin, $pos, 1);
    }
    return $string || "_";
}


sub help {
print STDERR "

This script computes a histogram of the distance between the 5' and 3' end of each strand to 3' and 5' end of the closest complementary strand. It also computes the closest distance between 5'/3' ends that are on the same strand for comparison. If desired, strands that are in distance X to the closest 3' or 5' end of the complementary strand can be written to separate bam files. Only alignments to the autosomes or chromosome X are analyzed.

CAUTION: For mate pairs, only the first read of the pair is currently written to bam output. Histograms are computed correctly. 

              5'                3'                    5'             3'                5'              3'     5'            3'
              ------Seq1-------->                 --------Seq2------->                 -----Seq3------->      -----Seq4----->
            <-----------------------                 <-----------        <--------                            <--------------
            3'                     5'                3'         5'       3'      5'                           3'            5'

5' distance:         -2                                    +3                                -14                     0
            <>                                    <->                    <------------>                       |

3' distance:         +3                                    -5                                +21                     0
                                 <->                             <--->                                  <------------------->

5'control            +36                                   +36                               +22                     +22
distance      <---------------------------------->                                      <-------------------->


[usage]
./complementary_strands2.pl [-option] in.bam

[options]
-quality      map quality cutoff [0]
-minread      minimal length cutoff [30]
-range        range of values to report [100]
-output       write out sequences that have another sequence in distance X from their 5' or 3' end respectively. 
-control      write out sequences that have another sequence in distance X from their 5' or 3' end respectively IN THE SAME ORIENTATION
-mt           analyze mtDNA (drop mapping restrictions to X and autosomes)

[output]
-in.hist.txt  table with histogram of closest distances, in number of observations and fraction of sequences (>redirect to file)
-in.ends.txt  table showing the observed combinations between 5' and 3' molecule structures
-in.matching_hist.txt table providing the fraction of sequences with a maximum distance of X between complementary strands, binned by sequence length 
-bam files    in.{5p/3p}-{overhang/recessed/blunt}-[X]bp.bam

";
exit;
}
