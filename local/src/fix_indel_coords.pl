#!/usr/bin/perl

use warnings;
use strict;

# write something relevant here TODO
use Getopt::Long;

my $help=0;
my $snp_info="";
my $len=0;

my $usage=<<END;
Usage: $0 -s snp_info < fimo output
	-s | --snp snp_info: indels info, format is rs, length, ins|del
    -l | --length: length of sequence before mutations
END

GetOptions (
	'help'=>\$help,
	'snp=s' =>\$snp_info,
	'length=i' =>\$len,
) or die($usage);


if($help){
	print($usage);
	exit(0);
}

if($snp_info eq "" || $len==0){
	print "A mandatory argument is missing!\n";
	die($usage);
}

open(SNP, $snp_info) || die "Error: $!, opening $snp_info\n";

my %indels;
while (<SNP>) {
    chomp;
    my @line = split("\t", $_);
    if (defined $indels{$line[0]}) {
        die "Error duplicated indel $line[0]\n";
    } else {
        $indels{$line[0]} = [$line[1], $line[2]]; #line[2] is ins or del, line[1] is the length  
    }
}

while(<STDIN>) {
    #[egrassi@gncomp3 huang_inhouse_motifbr]$ head lead/fimo_indels.txt.tmp 
    #FOXL1_HUMAN.H10MO.D             het     allele1 rs146029108     chr9    139329936       139330002       rs146029108-    28      46      +       20.7097 7.53e-08        0.00157TAGTTATTTATTTATTTTT
    chomp;
    my $coord;
    my @line = split("\t", $_);
    if ($line[3] eq "allele2") {
        my $rs = $line[4];
        $coord = $line[9];
        my $indel = $indels{$rs};
        my $indel_len = $indel->[0] -1;
        my $end_indel = $len+$indel_len;
        if ($indel->[1] eq "ins") {
            # we leave as they are coords <= $len, cause the indel does not disrupt them
            # fimo reports 1 based coords, not reversed for - strand matches.
            if ($coord > $len) {
                # inside the INdels we use fake coords to avoid matching the matches (sic) starting here on the alt allele.
                if ($coord <= $end_indel) {
                    $coord = "I" . ($end_indel-$coord);
                } else {
                    # after the INdels we need to change coords to match right matches (even if they will never be motif changing right now)
                    $coord = $coord - $indel_len;
                } 
            }
        } elsif ($indel->[1] eq "del") {
            # we leave as they are coords <= $len, cause the indel does not disrupt them
            # fimo reports 1 based coords, not reversed for - strand matches.
            if ($coord > $len) {
                # after the inDELS we need to change coords to match right matches (even if they will never be motif changing right now)
                $coord = $coord + $indel_len;
                } 
        } else {
           die "Error for $rs in info! $indel->[1]\n"
        }
        print $line[0] . "\t\t" . $line[2] . ';' . $line[3] . ';' . $line[4] . ';' . $line[5] . ';' . $line[6] . ';'  . $line[7] . ';' . $line[8] . "\t";
        print $coord . "\t" . $line[10] . "\t" . $line[11] . "\t" . $line[12] . "\t" . $line[13] . "\t" . $line[14] . "\t" . $line[15] . "\n";
    } else {
        print $line[0] . "\t\t" . $line[2] . ';' . $line[3] . ';' . $line[4] . ';' . $line[5] . ';' . $line[6] . ';'  . $line[7] . ';' . $line[8] . "\t";
        print $line[9] . "\t" . $line[10] . "\t" . $line[11] . "\t" . $line[12] . "\t" . $line[13] . "\t" . $line[14] . "\t" . $line[15] . "\n";
    }
}
