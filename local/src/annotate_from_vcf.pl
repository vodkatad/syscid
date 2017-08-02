#!/usr/bin/perl

use warnings;
use strict;

use lib '/data/egrassi/snaketree/prj/syscid/local/src/vcftools-vcftools-1d27c24/src/perl';
use Getopt::Long;
use Vcf;

my $help=0;
my $snps="";
my $info="";

my $usage=<<END;
Usage: $0 -s snps -i info < vcf
	-s | --snp: rs ids
    -i | --info: wanted info

END

GetOptions (
	'help'=>\$help,
	'snp=s' =>\$snps,
	'info=s' =>\$info,
) or die($usage);


if($help){
	print($usage);
	exit(0);
}

if($snps eq "" || $info eq ""){
	print "A mandatory argument is missing!\n";
	die($usage);
}

open(SNP, $snps) || die "Error: $!, opening $snps\n";
my %rs_ids;
my $locus_id = 0;
while (<SNP>) {
    chomp;
    my @line = split("\t", $_);
    my @rss = split("_", $line[3]);
    foreach my $rs (@rss) {
        if (defined $rs_ids{$rs}) {
            die "Error duplicated snp $_\n";
        } else {
            $rs_ids{$rs} = "$locus_id\t$line[0]\t$line[1]\t$line[2]";
        }
    }
    $locus_id++;
}
close(SNP);

open(INFO, $info) || die "Error: $!, opening $info\n";
my @wanted_info;
while (<INFO>) {
    chomp;
    push(@wanted_info, $_);
}
close(INFO);

my $vcf = Vcf->new(fh=>\*STDIN);
$vcf->parse_header();
my $ordered_id = 1;
while (my $x=$vcf->next_data_hash()) {
	my $id = $x->{"ID"};
    my $locus_info = $rs_ids{$id};
    if (defined($locus_info)) {
    	my $chr = $x->{"CHROM"};
    	my $coord = $x->{"POS"};
    	my $reference = $x->{"REF"};
	    my $alts = $x->{"ALT"};
        my $alt = "$alts->[0]";
        for (my $i = 1; $i < scalar( @{ $alts } ); $i++) {
            $alt = $alt . "," . $alts->[$i];
        }
        print "$locus_info\t$chr\t$coord\t$id\t$reference\t$alt";
        foreach my $i (@wanted_info) {
            my $value = $x->{"INFO"}{$i};
        	if (defined($value)) {
                print "\t$value"; 
            } else {
                print "\tNA";
            }
        }
        print "\n";
    }
}
