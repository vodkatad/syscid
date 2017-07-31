#!/usr/bin/perl

use warnings;
use strict;

# Given a list of SNPs and a file with LD info (haploreg format) produces a format suitable to get loci with ld_merge.py, i.e
# it adds coordinates and puts all pairs of snps on a single line.
# output: chr coord_1 coord_2 snp_1 snp_2 r2
use Getopt::Long;

my $help=0;
my $snps="";
my $coords="";

my $usage=<<END;
Usage: $0 -s snps -c cutoff < ld_info_haploreg_format (rs10    rs2237570,0.3,-0.72; ...)
	-s | --snp: rs ids
    -c | --coords: coords for SNPs
    output: chr coord_1 coord_2 snp_1 snp_2 r2
END

GetOptions (
	'help'=>\$help,
	'snp=s' =>\$snps,
	'coords=s' =>\$coords, # 1,2,3 from vcf
) or die($usage);


if($help){
	print($usage);
	exit(0);
}

if($snps eq "" || $coords eq ""){
	print "A mandatory argument is missing!\n";
	die($usage);
}

open(SNPS, $snps) || die "Error: $!, opening $snps\n";

my %rs_ids;
while (<SNPS>) {
    chomp;
    if (defined $rs_ids{$_}) {
        die "Error duplicated snp $_\n";
    } else {
        $rs_ids{$_} = 0;
    }
}
close(SNPS);

my %infos;
# note that haploregV4.0 vcf has hg19 coords (GRCh37.p5).
open(COORDS, $coords) || die "Error: $!, opening $coords\n";
while(<COORDS>) {
    chomp();
    my @line = split("\t", $_);
    my $chr = $line[0];
    my $coord = $line[1];
    my $rs = $line[2];
    # no check on repetitions
    $infos{$rs} = $chr . "_" . $coord;#[$chr, $coord];
}
close(COORDS);

while(<STDIN>) {
    #[egrassi@gncomp3 data]$ zcat LD_EUR.tsv.gz  | head -n 1
    #rs10    rs2237570,0.3,-0.72;rs117465896,0.75,-0.88;rs1972508,0.3,0.87;rs4727281,0.61,0.93;
    chomp;
    my @line = split("\t", $_);
    if (exists($rs_ids{$line[0]})) {
        $rs_ids{$line[0]} = 1;
        my $has_friends = 0;
        my $mine_info = $infos{$line[0]};
        if (!defined($mine_info)) {
            print STDERR "$line[0]\t" . 'nocoord' . "\n";
            next;
        } else {
            my @mine_infos = split('_', $mine_info);
            my @friends = split(';', $line[1]);
            foreach my $f (@friends) {
                my @r = split(',', $f);
                my $rs = $r[0];
                my $r2 = $r[1];
                my $friend_info = $infos{$rs};
                if (!defined($friend_info)) {
                    print STDERR "$rs\t" . 'fnocoord' . "\n";
                    next;
                } else {
                    $has_friends++;
                    my @friend_infos = split('_', $friend_info);
                    if ($friend_infos[0] ne $mine_infos[0]) {
                        die "Error for $line[0], $rs, $friend_info, $mine_info! Different chr?\n";
                    }
                    print "$mine_infos[0]\t$mine_infos[1]\t$friend_infos[1]\t$line[0]\t$rs\t$r2\n";
                }
            }
            if ($has_friends == 0) {
                # dummy printing
                print "$mine_infos[0]\t$mine_infos[1]\t0\t$line[0]\tdummy\t0\n";
            }
        }
    }
}

foreach my $rs (keys(%rs_ids)) {
    if ($rs_ids{$rs} == 0) {
        my $mine_info = $infos{$rs};
        if (!defined($mine_info)) {
            print STDERR "$rs\t" . 'nocoord' . "\n";
            next;
        } else {
            my @mine_infos = split('_', $mine_info);
            print "$mine_infos[0]\t$mine_infos[1]\t0\t$rs\tdummy\t0\n";
        } 
    }
}
