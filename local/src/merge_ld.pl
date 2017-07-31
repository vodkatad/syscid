#!/usr/bin/perl

use warnings;
use strict;

# Given a list of loci in the format locus_id rs_is merges loci with some SNPs in common.
use Getopt::Long;

my $help=0;

my $usage=<<END;
Usage: $0 < loci_by_find_ld
    output is locus_id \t rsid
END

GetOptions (
	'help'=>\$help,
) or die($usage);


if($help){
	print($usage);
	exit(0);
}

my %loci = {};
my %renamed_loci = {};
my $n_loci = 0;
while(<STDIN>) {
    chomp;
    my @line = split("\t", $_);
    my $rs = $line[1];
    my $locus = $line[0];
    if (exists($loci{$rs})) {
        # this locus needs a new name, but we cannot (or at least not easily) determine it now since we have not finished.
        push(@{ $rd_ids{$rs} }, $locus);
    }
    else {
        # cannot do the renaming at once otherwise we risk not renaming every locus (if we  
        #    if (exists($renamed_loci{$locus}) {
        #        $locus = $locus;
        #    }
        $rs_ids{$rs} = [$locus];
    }
}

foreach my $rs (keys(%loci)) {
    my @arr = @{ $loci{$rs} };
    if (scalar(@arr) >= 2) {
        my $new_name = join("_", @arr); # but later on what happens when we have a rs that again needs to rename one of this loci?
        # again it will be better to surrender an use graphs. We have a graph and we want the connected components.
    }
}
