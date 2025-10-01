#!/usr/bin/perl
# paste_fastq.pl
# by Nate Campbell
# paste R2 sequence to end of R1 in RC orientation for GT-seq primer analysis with PE data.
# Usage: provide R1.fastq and R2.fastq files as a command line arguments...

use strict; use warnings;

die "provide R1.fastq and R2.fastq files\n" unless @ARGV == 2;

my $file = "temp.fq";

`paste $ARGV[0] $ARGV[1] > $file`;

open (FASTQ, "<$file") or die "Error opening $file\n";

while (<FASTQ>) {
	my $infoline = $_;
	my $seqline = <FASTQ>;
	my $infoline2 = <FASTQ>;
	my $qual = <FASTQ>;
	my @infos = split "\t", $infoline;
	$infoline = $infos[0];
	chomp ($seqline);
	chomp ($qual);
	my @seqs = split "\t", $seqline;
	my @quals = split "\t", $qual;
	my $qual2 = reverse $quals[1];
	my $RC = reverse $seqs[1];
	$RC =~ tr/ACGT/TGCA/;
	$seqline = "$seqs[0]$RC";
	my $seqline2 = reverse $seqline;
	$seqline2 =~ tr/ACGT/TGCA/;
	$qual = "$quals[0]$qual2";
	my $qualRC = reverse $qual;
	print "$infoline\n$seqline\n+\n$qual\n";
	}
close FASTQ;

`rm $file`;
