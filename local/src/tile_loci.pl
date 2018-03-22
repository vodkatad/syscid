#!/usr/bin/perl

use warnings;
use strict;

# tile loci at a given tile lengths.
# we get the middle points for smaller intervals.
# we get a base every length*2 tiles otherwise (with the middle of the last smaller portion).
use Getopt::Long;

my $help=0;
my $tile=0;

my $usage=<<END;
Usage: $0 -l tile_length < merged_loci
	-l | --length: 
	TODO
END

GetOptions (
	'help'=>\$help,
	'length=i' =>\$tile,
) or die($usage);


if($help){
	print($usage);
	exit(0);
}

if($tile == 0){
	print "A mandatory argument is missing!\n";
	die($usage);
}
my $halftile = int($tile/2);

while(<STDIN>) {
	chomp;
	my @line = split("\t", $_);
	my $chr = $line[0];
 	my $b = $line[1];
	my $e = $line[2];
	my $snps = $line[3];
	my $snps_orig = $snps;
	my $length = $e - $b;
	if ($length <= $tile) {
		my $h = int(($e-$b)/2);
		&spit_out($chr, $b+$h, $snps);
	} else {
		my $ep = $b + $halftile;
		my $counter = 0;
		while ($ep < $e) {
			&spit_out($chr, $ep, $snps);
			$snps = $snps_orig . "-" . $counter;
			$counter++;
			$ep += $tile;
		}
		#my $covered = $counter*$tile; # -1 or not?
		my $reached = $ep-$tile;
		$ep = $reached+int(($e-$reached)/2);
		if ($ep > $e) {
			die "Something bad for $snps, $ep, $e, $tile";
		}
		&spit_out($chr, $ep, $snps);
	}
}

sub spit_out {
	my ($chr, $start, $id) = @_;
	#print "$chr. $start . $id \n";
	print 'chr'.$chr . "\t" . $start . "\t" . ($start+1) . "\t" . $id . "\n";
}
