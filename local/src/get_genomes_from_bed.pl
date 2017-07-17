#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use diagnostics;
$"="\t";

$SIG{__WARN__} = sub {die @_};

my $usage=<<END;
Usage: $0 [-c] [-i] sorted_alleles seq_dir bed_file < sorted_alleles_samples >
	-a|add_string		prepends the given string to the chr ids when looking for .fa in seq_dir
	-r|as_reference		converts dots (.) in the sorted_alleles_samples file as the reference genome variant. Otherwise an N will be inserted. Suggested option! 
	-o|overlapping_bed	doesn't check for overlapping bed, bed must be sorted. Default checks
	-e|empty_bed_file	checks if bed_file is empty and give a warning if this happens. Default gives no output and crashes
	
Will work on single chromosome basis.

bed files must be 0 based. Coordinates for SNP must be 1 based. Output is 0 based.

.META: sorted_alleles
	1	SNP_ID	snp_22_16057210
	2	chr	22
	3	pos	16057210
	4	ref	C
	5	alt	T

.META: stdin
	1	SNP_ID		snp_22_16057210
	2	SAMPLE_ID	HG00105
	3	allele1		T
	5	allele2		C

Example:
END
 
my $help=0;
my $ignore_missing_snp=0;
my $add_string="";
my $as_reference=0;
my $ignore_acrossbed_snp=0;
my $overlapping_bed=0;
my $empty_bed_file=0;
GetOptions (
	'h|help' => \$help,
	'a|add_string=s' => \$add_string,
	'r|as_reference' => \$as_reference,
	'o|overlapping_bed' => \$overlapping_bed,
	'e|empty_bed_file' => \$empty_bed_file
) or die($usage);

if ($help){
	print $usage;
	exit(0);
}

if ($#ARGV != 2) {
	die($usage);
}

my $alleles_info_file = $ARGV[0];
my $seq_dir = $ARGV[1];
my $bed_file = $ARGV[2];
my $reference_file = $add_string;

open(my $alleles_info_fh, "<", $alleles_info_file) || die "can't open $alleles_info_file: $!";

my @coords = ();
my @lengths_ref = ();
my @lengths_alt = ();
my @SNP_reference = ();
my @SNP_alternative = ();
my @snp_names = ();
#my $snp_id = 0;
my $coord = 1;
my $variant = "";
my $chr = "";
my $id = 0;
### OLD COMMENTS ###
# We have two problems, snp like these:
#snp_9_217892    1532    9       217892  G       .
#sv:4I_9_125223  640     9       125223  T       <DEL>
# The deletion following spec have info on length in the original vcf but they are lost during deconvolution. Just remove
# them if they are not too many?
# The dot are unknown values, so why they are there? In two examples they were present in the allele file but not
# in the samples, but who knows?
### OLD COMMENTS ###

# Here we read the allele info file
while(<$alleles_info_fh>) {
	chomp;
	my @line = split("\t",$_);
	# 2 is chr, 1 is ordered_snp_id, 3 is coords, 4 is reference seq (used to get length).
	# Order check: the file should be sorted on coords (corresponding to our internal numeric ids). 
	if ($id >= $line[1] || $coord > $line[3]) {
		die "Warning! Duplicate variant $id or alleles_info_file not sorted on coords\n";
	}
	if ($chr ne $line[2] && $chr ne "")  {
		die "Warning! More than one chr\n";
	}
	my $snp_id = $line[0];
	$id = $line[1];
	$chr = $line[2];
	$coord = $line[3];
	my $ref = $line[4];
	my $alt = $line[5];
	# We cannot push as long as we want to support situation in which some snp found here are not
	# in the samples file.
	# VCF files are 1 based so we convert coords to 0 based
	# The whole script should be 0 based. Please check this out!
	$coords[$id-1] = $coord-1;
	$lengths_ref[$id-1] = length($ref);
#	if ($alt =~ /DEL/) {
#               my @DEL = split("_", $alt);
#                $lengths_alt[$id-1] = $DEL[1];
#        } else {
	$lengths_alt[$id-1] = length($alt);
#        }
	$SNP_reference[$id-1] = $ref;
	$SNP_alternative[$id-1] = $alt;
	$snp_names[$id-1] = $snp_id;
	#push(@lengths_alt, length($alt));
}
$reference_file .= $chr;

close($alleles_info_fh) || warn "close failed: $!";

# Read reference genome in a whole string.
$reference_file .= ".fa";
open(my $reference_fh, "<", $seq_dir . "/"  . $reference_file) || die "Error opening $reference_file file: $!"; 
my $reference = "";
unless (<$reference_fh> =~ /^>/) {die "Reference file must be in fasta format"}
while(<$reference_fh>) {
#	next if /^>/; 
	# FIXME forzare ad avere solo prima riga con > e poi accettare solo ACTGN
	chomp;
	unless ($_ =~ /[ACTGN]/gi) {
		die "Reference sequence must be a fasta file and must contain only ACTGNactgn characters";
	}
	$reference .= $_;
}

close($reference_fh);

# Read .bed file as a whole
open(my $bed_fh, "<", $bed_file) || die "Error in opening $bed_file file: $!";
my @bed_coords = ();
my $bed_index = 0;
my $old_start_coord = -1;
my $old_end_coord = -1;
while(<$bed_fh>){
	chomp;
	my @line = split("\t", $_);
	# Check for unsorted regions
	if ($old_start_coord > $line[1]) {
		die "Warning! segments defined in bed file are not ordered. Bed file must be ordered!"
	}
	#Check for overlapping regions
	if (!$overlapping_bed && $old_end_coord >= $line[1]) {
		die "Warning! segments defined in bed file are overlapping. Segment: \"" . "@{$bed_coords[$bed_index-1]}" . "\"\toverlaps with: \"". "@line". "\"\n";
	}
	$bed_coords[$bed_index] = \@line;
	$old_start_coord = $line[1];
	$old_end_coord = $line[2];
	$bed_index++;
}

close($bed_fh);



# Now we loop over samples.
my $sample_id = "";
my $old_coord = 0;
# vcf have 1 based coords but we use 0 based coords throughout the script
my $current_ref_pos = 0;
my $bed_array_pos = 0;
my ($name_allele1, $name_allele2);
my ($bed_line, $next_bed_b);
my $old_sample_id = "";
my @snp_buffer = ();
my $old_snp_id = 0;
# We initialize bed variables for the first bed region
my $bed_chr = ${$bed_coords[0]}[0];
my $bed_b = ${$bed_coords[0]}[1];
my $bed_e = ${$bed_coords[0]}[2];
my $bed_name = ${$bed_coords[0]}[3];
# Checks if bed_file is empty and if so process ends here
my $is_ther_any_bed_region = 1;
if (!$bed_chr && !$bed_b && !$bed_e && !$bed_name && $empty_bed_file) {
	print STDERR "\nWARNING! bed_file: " . $bed_file . " is empty!\n";
	$is_ther_any_bed_region = 0;
}

if ($is_ther_any_bed_region) {
while(<STDIN> ) {
        chomp;
	#Read input SNP
# $line = [0] internal_id; [1] sample_id; [2] allele1; [3] allele2; [4] snp_coord; [5] length_ref; [6] length_alt; [7] snp_reference; [8] snp_names
        my @line = split("\t",$_);
	my $snp_id = $line[0];
	my $snp_coord = $coords[$snp_id-1];
#	my $length_ref = $lengths_ref[$snp_id-1];
#	my $length_alt = $lengths_alt[$snp_id-1];
	# FIXME CHECK ORDER degli snp_id e sample_id
	if (($old_sample_id eq $line[1]) && ($old_snp_id > $snp_id)) { die "SNP must be sorted on SNP ID! \n\n$!"}
	if ($old_sample_id gt $line[1]) { die "Samples must be sorted on standard input!\n\n$!"}
#	push (@line, $snp_coord);
#	push (@line, $length_ref);
#	push (@line, $length_alt);
#	push (@line, $SNP_reference[$snp_id-1]);
#	push (@line, $snp_names[$snp_id-1]);
#	push (@line, $SNP_alternative[$snp_id-1]);
	my ($length_ref, $length_alt) = 0;
	# Here we modify variants with DEL to be correctly processed downstream
	if ($SNP_alternative[$snp_id-1] =~ /DEL/) {
		my @DEL = split("_", $SNP_alternative[$snp_id-1]);
		$length_ref = $DEL[1];
		$length_alt = 1;
		for (my $allele_seq = 2; $allele_seq < 4; $allele_seq++) {
			if ($line[$allele_seq] =~ /DEL/) {
				$line[$allele_seq] = $SNP_reference[$snp_id-1];
#	        	        $line[6] = $DEL[1]+1; # XXX In caso di deletion avanziamo di $line[6], altrimenti di $line[5].
						      # Il DEL_X e` sempre nell'alternative.
			} else {
				$line[$allele_seq] = substr($reference, $snp_coord, $DEL[1]);
			}
	        }
	} else {
	        $length_ref = $lengths_ref[$snp_id-1];
	        $length_alt = $lengths_alt[$snp_id-1];
	}
	push (@line, $snp_coord);
	push (@line, $length_ref);
        push (@line, $length_alt);
        push (@line, $SNP_reference[$snp_id-1]);
        push (@line, $snp_names[$snp_id-1]);
	$old_snp_id = $snp_id;
#print STDERR "@line, ".$le."\n";
	$sample_id = $line[1];
	# Here we change sample, so snp_buffer and bed_array_pos must be set to 0
	if (($old_sample_id ne $sample_id) && ($old_sample_id ne "")) {
		# The last bed with SNP for the previous sample is printed here
#		&print_bed_sequence($old_sample_id, $bed_chr, $bed_name, $bed_b, $bed_e, \@snp_buffer, \$reference, $as_reference);
#		@snp_buffer = ();
#		my $snp_buffer = &update_snp_buffer($bed_chr,$bed_b,$bed_e,\@snp_buffer);
#               @snp_buffer = @{$snp_buffer};
		&print_bed_sequence($old_sample_id, $bed_chr, $bed_name, $bed_b, $bed_e, \@snp_buffer, \$reference, $as_reference);
		#print STDERR "ho finito con questo campione, ne inizio un altro! e passo qui otto volte\n";
		# Here we consider bed without SNP at the end of the chromosome and print them
		while ($bed_coords[$bed_array_pos+1]) { 
			$bed_array_pos++;
			$bed_line = $bed_coords[$bed_array_pos];
			$bed_chr = ${$bed_line}[0];
                        $bed_b = ${$bed_line}[1];
                        $bed_e = ${$bed_line}[2];
                        $bed_name = ${$bed_line}[3];
			my $snp_buffer = &update_snp_buffer($bed_chr,$bed_b,$bed_e,\@snp_buffer);
		        @snp_buffer = @{$snp_buffer};
			&print_bed_sequence($old_sample_id, $bed_chr, $bed_name, $bed_b, $bed_e, \@snp_buffer, \$reference, $as_reference);
		}
		# Resetting bed coordinates
		$bed_array_pos = 0;
		$bed_line = $bed_coords[$bed_array_pos];
		$bed_chr = ${$bed_line}[0];
                $bed_b = ${$bed_line}[1];
                $bed_e = ${$bed_line}[2];
                $bed_name = ${$bed_line}[3];
		@snp_buffer = ();
	}
	# Here we check if the current SNP overlaps with the previous. In this case we skip it
	# XXX TODO CHECK OVERLAP conditions: dobbiamo buttare via solo quelli che si overlappano sulla reference, gli
	# altri hanno senso (controllare i casi)
	if (($snp_buffer[-1]) && (((${$snp_buffer[-1]}[4] + ${$snp_buffer[-1]}[5]) > $line[4]) )) { # >= and not > is not the bug I am looking for (and > is ok)
		$old_sample_id = $sample_id;
		next;
	}
#	$old_snp_id = $snp_id;
	# If SNP is inside the bed insert it into SNP buffer
	if (($snp_coord + $length_ref > $bed_b) && ($snp_coord < $bed_e)) {
		push (@snp_buffer, \@line); 
		$old_sample_id = $sample_id;
#print STDERR "Actual snp buffer: @{$snp_buffer[0]}\n";
	# If SNP is outside the bed, print the bed and change bed
	} elsif ($snp_coord >= $bed_e) {
		# Here we print bed without SNP before the SNP we are now considering until the current SNP falls inside a bed region or before a bed region
		while ($snp_coord >= $bed_e && ($bed_coords[$bed_array_pos+1])) {
#print STDOUT "printing bed with snp!\n";
#			@snp_buffer = ();
			&print_bed_sequence($sample_id, $bed_chr, $bed_name, $bed_b, $bed_e, \@snp_buffer, \$reference, $as_reference);
			$bed_array_pos++;
			#next unless($bed_coords[$bed_array_pos]); # XXX toglibile?
			$bed_line = $bed_coords[$bed_array_pos];
			$bed_chr = ${$bed_line}[0];
			$bed_b = ${$bed_line}[1];
			$bed_e = ${$bed_line}[2];
			$bed_name = ${$bed_line}[3];
			my $snp_buffer = &update_snp_buffer($bed_chr,$bed_b,$bed_e,\@snp_buffer);
			@snp_buffer = @{$snp_buffer};
#print STDERR "updated snp buffer:  @$snp_buffer\n";
			# If SNP is inside new bed instert it into SNP buffer, else go to next SNP but retain this bed
			if (($snp_coord + $length_ref > $bed_b) && ($snp_coord < $bed_e)) {
				push (@snp_buffer, \@line);
			}
		}
		$old_sample_id = $sample_id;
	} else {
		# Skip SNP outside bed regions (in this case before current bed region)
		$old_sample_id = $sample_id;
		next;
	}
#print STDERR "Actual snp buffer: @{$snp_buffer[0]}\n";
}

# Here we print the last bed without SNP for the last sample
&print_bed_sequence($sample_id, $bed_chr, $bed_name, $bed_b, $bed_e, \@snp_buffer, \$reference, $as_reference);
#print STDERR "Actual snp buffer: @{$snp_buffer[0]}\n";
#@snp_buffer = ();
#my $snp_buffer = &update_snp_buffer($bed_chr,$bed_b,$bed_e,\@snp_buffer);
#@snp_buffer = @{$snp_buffer};
#print STDOUT "Ho finito col ciclo sullo STDIN!!\n";
while ($bed_coords[$bed_array_pos+1]) {
#print STDERR "Printing bed without snp\n";
	$bed_array_pos++;
	$bed_line = $bed_coords[$bed_array_pos];
	$bed_chr = ${$bed_line}[0];
	$bed_b = ${$bed_line}[1];
	$bed_e = ${$bed_line}[2];
	$bed_name = ${$bed_line}[3];
	my $snp_buffer = &update_snp_buffer($bed_chr,$bed_b,$bed_e,\@snp_buffer);
	@snp_buffer = @{$snp_buffer};
	&print_bed_sequence($sample_id, $bed_chr, $bed_name, $bed_b, $bed_e, \@snp_buffer, \$reference, $as_reference);
}
}

sub update_snp_buffer{
	my $bed_chr = shift;
        my $bed_b = shift;
        my $bed_e = shift;
        my $snp_buffer = shift;
	my @snp_buffer = @{$snp_buffer};
	my @new_snp_buffer = ();
	for (my $i = 0; $i < scalar(@snp_buffer); $i++) {
		my $snip_id = ${$snp_buffer[$i]}[0];
		my $sample_id = ${$snp_buffer[$i]}[1];
		my $snp_coord = ${$snp_buffer[$i]}[4];
		my $snp_ref_length = ${$snp_buffer[$i]}[5];
		if (($snp_coord + $snp_ref_length > $bed_b) && ($snp_coord < $bed_e)) {
			push (@new_snp_buffer, $snp_buffer[$i]);
		}
	}
	@snp_buffer = @new_snp_buffer;
#print STDERR "This is the updated snp_buffer: "." @$snp_buffer[0]";
	return(\@snp_buffer);
}
    
sub print_bed_sequence {
	my $sample_id = shift;
	my $bed_chr = shift;
	my $bed_name = shift;
	my $bed_b = shift;
	my $bed_e = shift;
	my $snp_buffer = shift;
	my @snp_buffer = @{$snp_buffer};
	my $reference = shift;
	my $as_reference = shift;
	my $allele = 2;
	my $all_snp_ids = "";
	if (@snp_buffer) {
		foreach my $s (@snp_buffer) {
			$all_snp_ids .= ${$s}[8] . "-";
		}
	}
	for (my $i = 1; $i < 3; $i++) {
		my $current_ref_pos = $bed_b;
		print STDOUT ">".$sample_id.";allele".$i.";".$bed_name.";".$bed_chr.";".$bed_b.";".$bed_e.";".$all_snp_ids."\n";
		if (@snp_buffer) {
#print STDERR $sample_id.";allele".$allele.";".$bed_name." contains SNP\n";
			foreach my $snp (@snp_buffer) {
				my $snp_id = ${$snp}[0];
		                my $sample_id = ${$snp}[1];
				my $allele_seq = "";
				$allele_seq = ${$snp}[$allele];
		                my $snp_coord = ${$snp}[4];
		                my $snp_ref_length = ${$snp}[5];
				my $snp_alt_length = ${$snp}[6];
				my $reference_seq = ${$snp}[7];
				# Treating special cases
			        if ($as_reference) {
			                if ($allele_seq =~ /^\.$/) {$allele_seq = $reference_seq}
			        } else {
			                if ($allele_seq =~ /^\.$/) {$allele_seq = "N"}
			        }
				if ($current_ref_pos < $snp_coord) {
					$current_ref_pos += &print_unchanged_seq($reference, $current_ref_pos, $snp_coord);
					$current_ref_pos += &print_snp ($snp_ref_length, $allele_seq, $sample_id, $reference_seq);
				} else {
					&print_snp ($snp_ref_length, $allele_seq, $sample_id, $reference_seq);
					$current_ref_pos = $snp_coord + $snp_ref_length;
				}
#				$current_ref_pos += &print_snp ($snp_ref_length, $allele_seq, $sample_id, $reference_seq);
			}
			if ($current_ref_pos < $bed_e) {
				&print_unchanged_seq($reference, $current_ref_pos, $bed_e);
			}
		} else {
#print STDERR $sample_id.";allele".$allele.";".$bed_name.";".$length_ref." don't contains SNP\n";
			&print_unchanged_seq($reference, $bed_b, $bed_e);
		}
		print STDOUT "\n";		
		$allele++;
	}
}

sub print_unchanged_seq {
        # mind example with a snp in coord 10.
        # _____________________
        # 1 2 3 4 5 6 7 8 9 | 10
        my $reference_seqR = shift;
        my $current_ref_pos = shift;
        my $snp_coord = shift;
        my $diff = $snp_coord - $current_ref_pos;
        if ($diff < 0) {
                die "WARNING! I am printing a \"negative\" substring! You probably encountered a special case not considered in this script. Debugging is absolutely required!";
        }
        my $to_print = substr($$reference_seqR, $current_ref_pos, $diff);
        print STDOUT $to_print;
        return($diff);
}

sub print_snp {
        #my $length_alt = shift;
        my $length_ref = shift;
        my $allele_seq = shift;
	my $sample_id = shift;
	my $reference_seq = shift;
#	if ($allele_seq =~ /DEL/) {
#		my @DEL = split("_", $allele_seq);
#		$length_ref = $DEL[1]+1;
#		$allele_seq = $reference_seq;
#	}
	print STDOUT $allele_seq;
        return($length_ref);
}

