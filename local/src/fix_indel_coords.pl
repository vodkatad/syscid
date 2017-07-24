#!/usr/bin/perl

use warnings;
use strict;

# write something relevant here TODO
use Getopt::Long;

my $help=0;
my $snp_info="";

my $usage=<<END;
Usage: $0 -s snp_info < fimo output
	-s | --snp snp_info: indels info, format is rs, length, ins|del
END

GetOptions (
	'help'=>\$help,
	'snp=s' =>\$snp_info,
) or die($usage);


if($help){
	print($usage);
	exit(0);
}

if($snp_info eq ""){
	print "A mandatory argument is missing!\n";
	die($usage);
}

# salvare in un array gli id che si desidera tenere
open (LISTA, $list) || die "Error: unable to open $list\n";
my @lista=<LISTA>;
chomp(@lista);
close LISTA;

# considerare la prima riga del file (header)
my $header=<STDIN>;
chomp($header);
my @header=split("\t",$header);

# stampare solo elementi selezionati e inseriti in @tenere
while(my $riga=<STDIN>){
	my @riga=split("\t",$riga);
	chomp(@riga);
	my @nuovariga=();
	for(my $i=0; $i<scalar(@tenere); $i++){
		push(@nuovariga, $riga[$tenere[$i]]);
	}
	# stampare la nuova riga
	print "@nuovariga\n";
} 
