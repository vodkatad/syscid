#!/usr/bin/perl

use warnings;
use strict;

#TODO
use Getopt::Long;
use Vcf;

my $help=0;
my $snps="";
my $cutoff=-1;

my $usage=<<END;
Usage: $0 -s snps -c cutoff < ld_info_haploreg_format (rs10    rs2237570,0.3,-0.72; ...)
	-s | --snp: rs ids
    -c | --cutoff: r2 threshold to consider two SNPS in LD

    output is locus_id \t rsid, the last one is the snp given in the snps file
END

GetOptions (
	'help'=>\$help,
	'snp=s' =>\$snps,
	'cutoff=f' =>\$cutoff,
) or die($usage);


exit(0);

if($help){
	print($usage);
	exit(0);
}

if($snps eq "" || $cutoff==-1){
	print "A mandatory argument is missing!\n";
	die($usage);
}

open(SNP, $snps) || die "Error: $!, opening $snps\n";

my %rs_ids;
while (<SNP>) {
    chomp;
    if (defined $rs_ids{$_}) {
        die "Error duplicated snp $_\n";
    } else {
        $rs_ids{$_} = 0;
    }
}

my $n_loci = 0;
while(<STDIN>) {
    #[egrassi@gncomp3 data]$ zcat LD_EUR.tsv.gz  | head -n 1
    #rs10    rs2237570,0.3,-0.72;rs117465896,0.75,-0.88;rs1972508,0.3,0.87;rs4727281,0.61,0.93;
    chomp;
    my @line = split("\t", $_);
    my $has_friends = 0;
    if (exists($rs_ids{$line[0]})) {
        my @friends = split(';', $line[1]);
        foreach my $f (@friends) {
            my @r = split(',', $f);
            my $rs = $r[0];
            my $r2 = $r[1];
            if ($r2 >= $cutoff) {
                $has_friends = 1;
                print "$n_loci\t$rs\n";
            }
        }
        #if ($has_friends) {
        print "$n_loci\t$line[0]\n";
        $n_loci++;
        #} #stupid Elena, y you want to forgot lonely singleton?
        $rs_ids{$line[0]} = 1;
    }
}

foreach my $rs (keys(%rs_ids)) {
    if ($rs_ids{$rs} == 0) {
        print "$n_loci\t$rs\n";
        $n_loci++;
    }
}
