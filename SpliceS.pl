#! /usr/bin/perl -w

use strict;
use Cwd;
use Cwd 'abs_path';
use IO::Handle;
use Getopt::Long;
use FileHandle;
use File::Basename;

##############################
# Arguments from command line
##############################

my $PROG = basename($0);
my $menu = Print_Menu();
my $command = '';
my $anno_rat = '';
my $anno_mouse = '';
my $anno_human = '';
my $input_file = '';
my $output_file = '';
my $known_species = 'rn,mm,hg';
my $novel_species = 'ps';

if(@ARGV == 0) { print "$menu\n"; exit(); }

$command = Get_Opt($command);

##############################
# Main body
##############################

print "\n[COMMANDS]\n$command\n";

print "\n[PROGRESS]\n";

my $seq_list;
my $anno_list;
my $anno_files = {
	'rn' => $anno_rat,
	'mm' => $anno_mouse,
	'hg' => $anno_human
};
my $exon_list;
my $splice_list;
my @known_species = split(/,/, $known_species);
my @novel_species = split(/,/, $novel_species);

Update_Progress('Parsing input file $input_file');
$seq_list = Parse_Fasta($input_file, $seq_list);

Update_Progress('Parsing gene annotation files');
foreach my $species (sort keys %{$anno_files}) {
	my $file = $anno_files->{$species};
	print "Reading [$species] $file ... \n";
	$anno_list = Parse_refFlat($file, $species, $anno_list);
}

# $anno_list = Print_Anno_Content($anno_list);

Update_Progress('Determining exon-exon boundries');
$exon_list = Determine_Exon($seq_list, $anno_list, $exon_list);

Update_Progress('Writing outputs');
$output_file = Write_Output($output_file);

Update_Progress("Program run finished\n");

##############################
# Functional Modules
##############################

sub Parse_Fasta 
{
	my ($input_file, $seq_list) = @_;

	my $func = 'Parse_Fasta';
	my $i = 0;
	my $line_skipped = 0;
	my $seq_total = 0;
	my $seq_header = '';
	my $seq = '';
	my $gene = '';
	my $acc = '';
	my $species = '';

	open(IN, '<', $input_file) or die $!;

	while(my $line = <IN>) {
		chomp($line);

		$line =~ s/\s+//g;

		$i++;

		if($i == 1 && $line !~ m/^>/) {
			print STDERR "$func: first line is not seq header. $input_file: $line. ",
			"Please check file format. Is this Fasta? ",
			"Program terminates. \n";
			exit(2);
		}

		if($line =~ m/^>/) {
			## this is the header [>gene|accessionID|species]
			$seq_total++;

			$seq_header = $line;
			$seq_header =~ s/^>//;

			my @header = split(/\|/, $seq_header);
			if(@header != 3) {
				print STDERR "$func: seq header not informat. $seq_header. ",
				"Program terminates. \n";
				exit(2);
			}

			($gene, $acc, $species) = @header;

			if($gene eq '' || $acc eq '' || $species eq '') {
				print STDERR "$func: gene, acc or species is empty. $seq_header. ",
				"Program terminates. \n";
				exit(2);
			}

			## format keys 
			$gene = lc($gene);
			$acc =~ s/\.\d+$//g;
			$species =~ s/\d+$//g; ## remove species genome versions (rn6 => rn)

			## add sequences
			if(exists $seq_list->{$gene}->{$acc}->{$species}) {
				print STDERR "$func: seq_list: key already exists.",
				"$gene=>$acc=>$species\n";
			} else {
				$seq_list->{$gene}->{$acc}->{$species} = '';
			}

		} elsif ($line =~ m/^[-ATGCN]+/) {
			## this is sequence line
			$seq_list->{$gene}->{$acc}->{$species} .= $line;

		} elsif ($line =~ m/\S+/) {
			print STDERR "$func: line not in format. $input_file: $line\n";
		}
	}

	close(IN);

	## print metrics 
	print "Line total = $i\nSeq total = $seq_total\n";
	# $seq_list = Print_Seq_Content($seq_list);

	return $seq_list;
}

sub Parse_refFlat
{
	my ($file, $species, $anno_list) = @_;

	my $func = 'Parse_refFlat';
	my $i = 0;
	my $line_skipped = 0;
	my $gene_total = 0;
	my $gene_recorded = 0;

	open(IN, '<', $file) or die $!;

	while(my $line = <IN>) {
		chomp($line);

		my @line = split(/\t/, $line);
		for(my $i = 0; $i < @line; $i++) { $line[$i] =~ s/\s+//g; }

		$i++;

		if($line =~ m/^\w+/) {
			$gene_total++;

			if(@line != 11) {
				print STDERR "$func: gene line not informat. $file: $line. ",
				"Please check file format. Is this refFlat? ",
				"Program terminates. \n";
				exit(2);
			}

			$line[0] = lc($line[0]);
			$line[1] =~ s/\.\d+$//g;

			my $gene = $line[0];
			my $acc = $line[1];

			$line = join("\t", @line);

			## no need to record all genes;
			## only record the gene if it is present in the input sequence file
			if(exists $seq_list->{$gene}->{$acc}->{$species}) {
				if(exists $anno_list->{$gene}->{$acc}->{$species}) {
					print STDERR "$func: anno_list: key already exists. ",
					"$gene=>$acc=>$species\n";
				} else {
					$anno_list->{$gene}->{$acc}->{$species} = $line;
					$gene_recorded++;
				}
			}			

		} elsif ($line =~ m/\S+/) {
			print STDERR "$func: line not in format. $input_file: $line\n";
		}
	}

	close(IN);

	## print metrics 
	print "Line total = $i\nGene total = $gene_total\n",
	"Gene recorded = $gene_recorded\n";

	return $anno_list;
}

sub Determine_Exon
{
	my ($seq_list, $anno_list, $exon_list) = @_;

	my $func = 'Determine_Exon';

	foreach my $gene (sort keys %{$seq_list}) {
		foreach my $acc (sort keys %{$seq_list->{$gene}}) {
			foreach my $species (sort keys %{$seq_list->{$gene}->{$acc}}) {
				my $seq = $seq_list->{$gene}->{$acc}->{$species};

				if($species eq 'ps') {
					next;
				}

				if(exists $anno_list->{$gene}->{$acc}->{$species}) {
					my @line = split(/\t/, $anno_list->{$gene}->{$acc}->{$species});

					$exon_list = Parse_Exon_Coordinate($exon_list, $species, \@line);

				} else {
					print STDERR "$func: anno_list: key does not exist. ",
					"$gene=>$acc=>$species\n";
				}
			}
		}
	}

	# $exon_list = Print_Exon_Content($exon_list);
	# $splice_list = Print_Splice_Content($splice_list);

	return $exon_list;
}

sub Write_Output
{
	my ($output_file) = @_;

	my $func = 'Write_Output';

	my $output_1 = "$output_file.seq.txt";
	my $output_2 = "$output_file.gene_model.txt";
	my $output_3 = "$output_file.splice_pos.txt";

	print "Printing $output_1 ...\n",
	"Printing $output_2 ...\n",
	"Printing $output_3 ...\n";

	open(OUT1, '>', $output_1) or die $!;
	open(OUT2, '>', $output_2) or die $!;
	open(OUT3, '>', $output_3) or die $!;

	my @header_out_2 = qw(GeneSymbol AccessionID Species 
	tx_start|tx_end exon_count exon_sizes intron_sizes cd_start|cd_end cd_sizes
	cd_exon_start_number|cd_exon_end_number);
	my @header_out_3 = qw(GeneSymbol AlnPos AccessionID Species 
	SplicePos IntronSize);

	print OUT2 join("\t", @header_out_2)."\n";
	print OUT3 join("\t", @header_out_3)."\n";

	foreach my $gene (sort keys %{$seq_list}) {
		## monitor if the current position is splice site in any sequence
		## if all pointers are 1, then this is conserved splice site 
		## across species
		my $pointer_list;
		my $splice_conserved; ## conserved splice sites

		## identify splice site position in aligned sequences
		$pointer_list = Scan_Sequence_Splice($gene, $seq_list, $pointer_list);

		## identify conserved splicing site
		# $pointer_list = Print_Pointer_Content($gene, $pointer_list);
		$splice_conserved = 
			Identify_Conserved_Splice($gene, $pointer_list, $splice_conserved);

		## print non-gapped sequence with splice position labeled
		$pointer_list = Write_Output_Sequence($gene, 
			$seq_list, $splice_list,$pointer_list);

		## print gene model info 
		$exon_list = Write_Output_Genemodel($gene, $exon_list);

		## print conserved splice info
		$pointer_list = Write_Output_Splice($gene, $pointer_list);
	}

	close(OUT1);
	close(OUT2);
	close(OUT3);

	return $output_file;
}

sub Identify_Conserved_Splice
{
	my ($gene, $pointer_list, $splice_conserved) = @_;

	my $func = 'Identify_Conserved_Splice';

	foreach my $pos (sort{$a<=>$b} keys %{$pointer_list}) {
		my $flag = 0;
		foreach my $species (@known_species) {
			if(exists $pointer_list->{$pos}->{$species} && 
				$pointer_list->{$pos}->{$species} =~ m/\d+/) {
				$flag++;
			}
		}

		if($flag == @known_species) {
			$splice_conserved->{$gene}->{$pos} = '';
		}
	}

	return $splice_conserved;
}

sub Scan_Sequence_Splice
{
	my ($gene, $seq_list, $pointer_list) = @_;

	my $func = 'Scan_Sequence_Splice';

	foreach my $acc (sort keys %{$seq_list->{$gene}}) {
		foreach my $species (sort keys %{$seq_list->{$gene}->{$acc}}) {
			my $seq = $seq_list->{$gene}->{$acc}->{$species};
			my $seq_nogap = $seq;
			$seq_nogap =~ s/-//g;
			# print "$gene, $acc, $species\n";

			my @chars = split(//, $seq);
			my $counter = 0;
			for(my $i = 0; $i < @chars; $i++) {
				my $char = $chars[$i];
				if($char =~ m/^[A-Z]$/i) {
					$counter++;
					
					if(exists 
						$splice_list->{$gene}->{$acc}->{$species}->{$counter}) {
						## does not consider multiple accs per gene!
						$pointer_list->{$i+1}->{$species} = $counter;
						# print substr($seq_nogap, 0, $counter)."\n";
					}

				} elsif ($char ne '-') {
					print STDERR "$func: char not in format. $char\n";
				}

			}
		}
	}

	return $pointer_list;
}

sub Parse_refFlat_Line
{
	my $line_ref = shift;

	my $func = 'Parse_refFlat_Line';
	my @line = @{$line_ref};
	# print "@line\n";
	if(@line != 11) {
		print STDERR "$func: this is not in refFlat format. line = \"@line\"\n";
	}

	my ($gene, $acc, $chr, $strand, 
		$tx_start, $tx_end, 
		$cd_start, $cd_end, $exon_count, 
		$exon_start, $exon_end) = @line;
	my @exon_starts;
	my @exon_ends;
	my @refflat_info;

	$exon_start =~ s/,$//;
	$exon_end =~ s/,$//;

	## refFlat has end coordinates open (1-offset, left close, right open)
	## minus all end position by 1
	$tx_end -= 1;
	$cd_end -= 1;
	@exon_ends = split(',', $exon_end);
	for(my $i = 0; $i < @exon_ends; $i++) { $exon_ends[$i] -= 1; }
	$exon_end = join(',', @exon_ends);

	@refflat_info = ($gene, $acc, $chr, $strand, $tx_start, $tx_end, 
		$cd_start, $cd_end, $exon_count,
		"$exon_start;$exon_end");

	return @refflat_info;
}

sub Parse_Exon_Coordinate
{
	my ($exon_list, $species, $line_ref) = @_;

	my $func = 'Parse_Exon_Coordinate';
	my ($gene, $acc, $chr, $strand, $tx_start, $tx_end,
		$cd_start, $cd_end, $exon_count, $exon) = Parse_refFlat_Line($line_ref);
	my @exons = split(/;/, $exon);
	my @exon_starts = split(/,/, $exons[0]);
	my @exon_ends = split(/,/, $exons[1]);
	my ($exon_starts_ref, $exon_ends_ref);

	# print "$gene | $acc | $chr | $strand | $tx_start | $tx_end | ",
	# 	"$cd_start | $cd_end | $exon_count | @exon_starts | @exon_ends\n";

	if(@exon_starts != @exon_ends || 
		@exon_starts != $exon_count ||
		@exon_ends != $exon_count) {
		print STDERR "$func: exon start/end is not consistent. ",
		"exon_count = $exon_count; exon_starts = $exons[0]; exon_ends = $exons[1]. ",
		"Program terminates.\n";
		exit(2);
	}

	## convert to relative genomic coordinates 
	if($strand eq '+' || $strand eq '-') {

		($tx_start, $tx_end, $cd_start, $cd_end, 
			$exon_count, $exon_starts_ref, $exon_ends_ref) = 
				Convert_Relative_Coordinate($exon_list, $species,
					$gene, $acc, $chr, $strand, $tx_start, $tx_end,
					$cd_start, $cd_end, $exon_count, 
					\@exon_starts, \@exon_ends);
		@exon_starts = @$exon_starts_ref;
		@exon_ends = @$exon_ends_ref;

		# print "$gene | $acc | $chr | $strand | $tx_start | $tx_end | ",
		# 	"$cd_start | $cd_end | $exon_count | @exon_starts | @exon_ends\n";

		$exon_list = Parse_Exon_Coordinate_Detail($exon_list, $species,
			$gene, $acc, $chr, $strand, $tx_start, $tx_end,
			$cd_start, $cd_end, $exon_count, 
			\@exon_starts, \@exon_ends);

	} else {
		print STDERR "$func: strand not in format. $strand\n";
	}

	return $exon_list;
}

sub Parse_Exon_Coordinate_Detail
{
	my ($exon_list, $species,
		$gene, $acc, $chr, $strand, $tx_start, $tx_end,
		$cd_start, $cd_end, $exon_count, 
		$exon_starts_ref, $exon_ends_ref) = @_;

	my $func = 'Parse_Exon_Coordinate_Detail';
	my @exon_starts = @$exon_starts_ref;
	my @exon_ends = @$exon_ends_ref;
	my @tx_locals = ('NA') x 2;
	my @cd_locals = (1, 'NA');
	my @cd_sizes = ();
	my @cd_exon_numbers = ('NA') x 2; 
	my @exon_sizes = ('NA') x $exon_count;
	my @intron_sizes = ('NA') x ($exon_count - 1);
	my $splice_pos = 0; 

	# print "$gene | $acc | $chr | $strand | $tx_start | $tx_end | ",
	# 	"$cd_start | $cd_end | $exon_count | @exon_starts | @exon_ends\n";

	## determine each exon length
	for(my $i = 0; $i < $exon_count; $i++) {
		$exon_sizes[$i] = $exon_ends[$i] - $exon_starts[$i] + 1;
	}

	## determine each intron length
	if($exon_count == 1) { $intron_sizes[0] = 0; } else {
		for(my $i = 1; $i < $exon_count; $i++) {
			$intron_sizes[$i-1] = $exon_starts[$i] - $exon_ends[$i-1] - 1;
		}
	}

	## determine whole transcript local start and end coordinates
	@tx_locals = (1, Calculate_Sum(\@exon_sizes));

	## determine CD local start and end on the transcript
	## first decide which exon has CD start

	for(my $i = 0; $i < $exon_count; $i++) {
		my $start = $exon_starts[$i];
		my $end = $exon_ends[$i];

		if($cd_start > $start && $cd_start > $end) {
			## have not reached CD start .... 
			$cd_locals[0] += $end - $start;

		} elsif($cd_start >= $start && $cd_end <= $end) {
			## the entire CD region is within one exon 
			$cd_locals[0] += $cd_start - $start + 1;
			push(@cd_sizes, $cd_end - $cd_start + 1);
			@cd_exon_numbers = ($i+1, $i+1);

		} elsif($cd_start >= $start && $cd_start <= $end && $cd_end > $end) {
			## reached the 1st CD exon, but not the last CD exon
			$cd_locals[0] += $cd_start - $start;
			push(@cd_sizes, $end - $cd_start + 1);
			$cd_exon_numbers[0] = $i+1;
		} elsif($cd_start < $start && $cd_end > $end) {
			## in the middle of CD exons
			push(@cd_sizes, $end - $start + 1);
		} elsif($cd_start < $start && $cd_end >= $start && $cd_end <= $end) {
			## reached the last CD exon
			push(@cd_sizes, $cd_end - $start + 1);
			$cd_exon_numbers[1] = $i+1;
		} else {
			print STDERR "$func: CD coordinates not considered. ",
			"exon number = $i\n";
		}
		
	}
	$cd_locals[1] = $cd_locals[0] + Calculate_Sum(\@cd_sizes) - 1;

	# print "@tx_locals; @cd_locals; @exon_sizes; @cd_sizes; @intron_sizes\n";
	# print "$gene, $acc, $species: [$strand] CD starts at exon $cd_exon_numbers[0] ",
	# "and ends at exon $cd_exon_numbers[1]\n";

	if(exists $exon_list->{$gene}->{$acc}->{$species}) {
		print STDERR "$func: exon_list: key already exists. ",
		"$gene=>$acc=>$species\n";
	} else {
		$exon_list->{$gene}->{$acc}->{$species} = join(';',
			join(',', @tx_locals), 
			$exon_count,
			join(',', @exon_sizes), 
			join(',', @intron_sizes),
			join(',', @cd_locals), 
			join(',', @cd_sizes), 
			join(',', @cd_exon_numbers)
			);
	}

	## record splice position in the sequences
	for(my $i = 0; $i < $exon_count - 1; $i++) {
		$splice_pos += $exon_sizes[$i];

		if(exists $splice_list->{$gene}->{$acc}->{$species}->{$splice_pos}) {
			print STDERR "$func: splice_list: key already exists. ",
			"$gene=>$acc=>$species=>$splice_pos\n";
		} else {
			$splice_list->{$gene}->{$acc}->{$species}->{$splice_pos} = 
				$intron_sizes[$i];
		}
		
	}
	
	return $exon_list;
}

sub Convert_Relative_Coordinate
{
	my ($exon_list, $species,
		$gene, $acc, $chr, $strand, $tx_start, $tx_end,
		$cd_start, $cd_end, $exon_count, 
		$exon_starts_ref, $exon_ends_ref) = @_;

	my $func = 'Convert_Relative_Coordinate';
	my @exon_starts = @$exon_starts_ref;
	my @exon_ends = @$exon_ends_ref;

	my $offset = '';
	if($strand eq '+') {
		$offset = $tx_start;

		$tx_start -= $offset + 1;
		$tx_end -= $offset + 1;
		$cd_start -= $offset + 1;
		$cd_end -= $offset + 1;
		for(my $i = 0; $i < $exon_count; $i++) {
			$exon_starts[$i] -= $offset + 1;
			$exon_ends[$i] -= $offset + 1;
		}
	} elsif($strand eq '-') {
		$offset = $tx_end;

		$tx_start = $offset - $tx_start + 1;
		$tx_end = $offset - $tx_end + 1;
		$cd_start = $offset - $cd_start + 1;
		$cd_end = $offset - $cd_end + 1;
		for(my $i = 0; $i < $exon_count; $i++) {
			$exon_starts[$i] = $offset - $exon_starts[$i] + 1;
			$exon_ends[$i]  = $offset - $exon_ends[$i] + 1;
		}

		## flip genomic coordinates for [-] strand
		($tx_start, $tx_end) = ($tx_end, $tx_start);
		($cd_start, $cd_end) = ($cd_end, $cd_start);
		my @tmp = @exon_starts;
		@exon_starts = reverse(@exon_ends);
		@exon_ends = reverse(@tmp);
	}

	return $tx_start, $tx_end, $cd_start, $cd_end, $exon_count, 
	\@exon_starts, \@exon_ends;
}

sub Write_Output_Sequence
{
	my ($gene, $seq_list, $splice_list,$pointer_list) = @_;

	my $func = 'Write_Output_Sequence';

	foreach my $acc (sort keys %{$seq_list->{$gene}}) {
		foreach my $species (sort keys %{$seq_list->{$gene}->{$acc}}) {
			my $seq_header = "$gene|$acc|$species";
			my $seq_out = '';
			my $seq = $seq_list->{$gene}->{$acc}->{$species};
			my $seq_nogap = $seq;
			$seq_nogap =~ s/-//g;
			# print "$gene, $acc, $species\n";

			my @chars = split(//, $seq);
			my $counter = 0;
			for(my $i = 0; $i < @chars; $i++) {
				my $char = $chars[$i];
				if($char =~ m/^[A-Z]$/i) {
					$counter++;
					$seq_out .= $char;

					## if this position is recorded as splice site
					## in the known species sequence
					if(exists $pointer_list->{$i+1}) {
						my @strings;
						foreach my $species (sort keys %{$pointer_list->{$i+1}}) {
							my $splice_pos = $pointer_list->{$i+1}->{$species};
							my $intron_size = 'NA';
							foreach my $acc (sort keys %{$splice_list->{$gene}}) {
								if(exists 
					$splice_list->{$gene}->{$acc}->{$species}->{$splice_pos}) {
								$intron_size = 
					$splice_list->{$gene}->{$acc}->{$species}->{$splice_pos};
								} 
							}
							
							push(@strings, "$species:$intron_size");
						}

						# print "@strings\n";
						$seq_out .= '['.join(',', @strings).']';
					}


				} elsif ($char ne '-') {
					print STDERR "$func: char not in format. $char\n";
				}
			}

			print OUT1 ">$seq_header\n$seq_out\n";
		}
	}

	return $pointer_list;
}

sub Write_Output_Genemodel
{
	my ($gene, $exon_list) = @_;

	foreach my $acc (sort keys %{$exon_list->{$gene}}) {
		foreach my $species (sort keys %{$exon_list->{$gene}->{$acc}}) {
			my $gene_model = $exon_list->{$gene}->{$acc}->{$species};
			$gene_model =~ s/;/\t/g;
			$gene_model =~ s/,/ \| /g;
			print OUT2 "$gene\t$acc\t$species\t$gene_model\n";
		}
	}

	return $exon_list;
}

sub Write_Output_Splice
{
	my ($gene, $pointer_list) = @_;

	foreach my $pos (sort{$a<=>$b} keys %{$pointer_list}) {
		foreach my $species (sort keys %{$pointer_list->{$pos}}) {
			my $splice_pos = $pointer_list->{$pos}->{$species};
			my $acc_id = 'NA';
			my $intron_size = 'NA';
			foreach my $acc (sort keys %{$splice_list->{$gene}}) {
				if(exists 
				$splice_list->{$gene}->{$acc}->{$species}->{$splice_pos}) {
					$acc_id = $acc;
					$intron_size = 
					$splice_list->{$gene}->{$acc}->{$species}->{$splice_pos};
				} 
				
			}

			print OUT3 "$gene\t$pos\t",
			"$acc_id\t$species\t$splice_pos\t$intron_size\n";
		}
	}

	return $pointer_list;
}

sub Print_Pointer_Content
{
	my ($gene, $pointer_list) = @_;

	# print "GeneSymbol\tAlnPos\tSpecies\tSplicePos\n";

	foreach my $pos (sort{$a<=>$b} keys %{$pointer_list}) {
		foreach my $species (sort keys %{$pointer_list->{$pos}}) {
			print "$gene | $pos | $species | ",
			$pointer_list->{$pos}->{$species}."\n";
		}
	}

	return $pointer_list;
}

sub Print_Exon_Content
{
	my $exon_list = shift;

	print "GeneSymbol\tAccessionID\tSpecies\t",
	"[tx_start,tx_end;exon_count;exon_sizes;",
	"intron_sizes;cd_start,cd_end;cd_sizes;",
	"cd_exon_start_number,cd_exon_end_number]\n";

	foreach my $gene (sort keys %{$exon_list}) {
		foreach my $acc (sort keys %{$exon_list->{$gene}}) {
			foreach my $species (sort keys %{$exon_list->{$gene}->{$acc}}) {
				my $line = $exon_list->{$gene}->{$acc}->{$species};
				print "$gene\t$acc\t$species\t[$line]\n";
			}
		}
	}

	return $exon_list;
}

sub Print_Splice_Content
{
	my $splice_list = shift;

	print "GeneSymbol\tAccessionID\tSpecies\tSplicePosition\tIntronSize\n";

	foreach my $gene (sort keys %{$splice_list}) {
		foreach my $acc (sort keys %{$splice_list->{$gene}}) {
			foreach my $species (sort keys %{$splice_list->{$gene}->{$acc}}) {
				foreach my $splice_pos (sort keys 
					%{$splice_list->{$gene}->{$acc}->{$species}}) {
					my $intron_size = 
					$splice_list->{$gene}->{$acc}->{$species}->{$splice_pos};
					print "$gene\t$acc\t$species\t$splice_pos\t$intron_size\n";
				}
			}
		}
	}

	return $splice_list;
}

sub Print_Seq_Content
{
	my $seq_list = shift;

	print "GeneSymbol\tAccessionID\tSpecies\tSeqLength\tSeqLengthNoGap\n";
	foreach my $gene (sort keys %{$seq_list}) {
		foreach my $acc (sort keys %{$seq_list->{$gene}}) {
			foreach my $species (sort keys %{$seq_list->{$gene}->{$acc}}) {
				my $seq = $seq_list->{$gene}->{$acc}->{$species};
				my $seq_nogap = $seq;
				$seq_nogap =~ s/-//g;
				print "$gene\t$acc\t$species\t".
					length($seq)."\t".length($seq_nogap)."\n";
			}
		}
	}

	return $seq_list;
}

sub Print_Anno_Content
{
	my $anno_list = shift;

	print "GeneSymbol\tAccessionID\tSpecies\trefFlat\n";
	foreach my $gene (sort keys %{$anno_list}) {
		foreach my $acc (sort keys %{$anno_list->{$gene}}) {
			foreach my $species (sort keys %{$anno_list->{$gene}->{$acc}}) {
				my $line = $anno_list->{$gene}->{$acc}->{$species};
				$line =~ s/\s+/ \| /g;
				print "$gene\t$acc\t$species\t\"$line\"\n";
			}
		}
	}

	return $anno_list;
}

sub Calculate_Sum
{
	my($data) = @_;

	my $func = 'Calculate_Sum';
	my $sum = 0;

	if (@$data == 0) {
	    print STDERR "$func: data array is empty.\n";
	} else {
		foreach (@$data) { $sum += $_ ; }
	}

	return $sum;
}

sub Calculate_Mean
{
	my($data) = @_;

	my $func = 'Calculate_Mean';
	my $sum = 0;
	my $mean = 'NA';

	if (@$data == 0) {
	    print STDERR "$func: data array is empty.\n";
	} else {
		foreach (@$data) { $sum += $_ ; }
		$mean = $sum / @$data;
	}

	return $mean;
}

sub Calculate_Stdev
{
	my($data) = @_;

	my $func = 'Calculate_Stdev';
	my $std = 'NA';
	my $mean = Calculate_Mean($data);
	my $sqtotal = 0;

	if (@$data == 0) {
	    print STDERR "$func: data array is empty.\n";
	} elsif (@$data == 1){
	    $std = 0;
	} else {
		foreach(@$data) {
			$sqtotal += ($mean - $_) ** 2;
		}
		$std = ($sqtotal / (@$data-1)) ** 0.5;
	}

	return $std;
}

sub Get_Opt
{
	my ($command) = @_;

	GetOptions( 
		'input|i:s' => \$input_file,
		'output|o:s' => \$output_file,
		'rat|rn:s' => \$anno_rat,
		'mouse|mm:s' => \$anno_mouse,
		'human|hg:s' => \$anno_human
	) or die $!; 

	if($input_file eq '' || 
		$anno_rat eq '' || $anno_mouse eq '' || $anno_human eq '') {
		print STDERR 'Input missing!\n$menu\n';
		exit(1);
	}

	if($output_file eq '') { $output_file = "$input_file.out"; }

	$command = "$PROG --input $input_file --output $output_file".
	"--rat $anno_rat --mouse $anno_mouse --human $anno_human";

	return $command;
}

sub Print_Menu
{
	my $sep = "-" x 80;
	my $version = '0.1';
	my $org = 'Center for Research Informatics, University of Chicago';

	my $menu = qq~
$sep
Copyright (c) 2016 $org

Usage: $PROG [OPTIONS]
Vesion: $version

Options:
 [-i|--input]     : Input sequence file. Fasta format.
                    Sequence header format must contain three fields,
                    separated by [|].  No space allowed. 
                    For example: >gene_symbol|transcript_accessionID|species
                    [species] must be 
                    rn or rn6; mm or mm10; hg or hg38.
             Note : All three homolog sequences are required (rn, mm, hg). 
 [-o|--output]    : Output file prefix. Default is "\$input.out".
 [-rn|--rat]      : Rat RefSeq gene annotation. refFlat format.
 [-mm|--mouse]    : Mouse RefSeq gene annotation. refFlat format.
 [-hg|--human]    : Human RefSeq gene annotation. refFlat format.
             Note : All three annotation files are required (rn, mm, hg). 

Three output files will be generated:
 \$output.seq.txt : Non-gapped transcript sequence with splice site labeled.
 \$output.gene_model.txt: Gene model information of known species.
 \$output.splice_pos.txt: Splice position information of the gene homologs.

Example: $PROG \\
	 -i KiSS1.hamster_RNAseq.alignment.fasta \\
	 -rn ../data/rn6.refFlat.txt \\
	 -mm ../data/mm10.refFlat.txt \\
	 -hg ../data/hg38.refFlat.txt

Contact: Riyue Bao <rbao\@uchicago.edu>
$sep
~;

	return $menu;
}

sub Parse_Inputfiles
{
	my @files = @_;

	my @files_new;
	
	if(@files == 1 && $files[0] =~ m/^(\S*)\*(\S*)$/){
		@files_new = <$1*$2>;
	} elsif (@files > 1){
		foreach my $file(@files) {
			if ($file =~ m/^(\S*)\*(\S*)$/) {
				@files_new = (@files_new, <$1*$2>);
			} else {
				push(@files_new, $file);
			}
		}
	} else {
		@files_new = @files;
	}

	return @files_new;
}


sub Update_Progress
{
	my $progress = shift;
	
	print "[",Local_Time(), "] $progress\n";
	
	return $progress;
}

sub Test_Hash
{
	my $hash = shift;
	
	foreach my $key (sort keys %{$hash}) {
		print "$key=>\n";
		
		if ($hash->{$key} =~ m/^HASH/) {
			Test_Hash($hash->{$key});
			print "\n";
		} elsif ($hash->{$key} =~ m/^ARRAY/){
			print join(" ", @{$hash->{$key}}), "\n";
		} else {
			print $hash->{$key}, "\n";
		}
	
	}

	return $hash;
}

##--- Print local time
sub Local_Time
{
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my  @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, 
		$yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	#my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek], $months[$month] $dayOfMonth, $year";
	
	my $digits = 2;
	$month = Add_Prefix($month, $digits);
	$dayOfMonth = Add_Prefix($dayOfMonth, $digits);
	$hour = Add_Prefix($hour, $digits);
	$minute = Add_Prefix($minute, $digits);
	$second = Add_Prefix($second, $digits);	
	my $theTime = "$month-$dayOfMonth-$year\_$hour:$minute:$second";
	
	#print "Local time: ", $theTime, "\n";
	return $theTime;
} 

sub Add_Prefix
{
	my ($number, $digits) = @_;
	
	if($number =~ m/^\d{$digits}$/){	
		return $number;
	}else{
		$number = "0".$number;
		Add_Prefix($number, $digits);
	}

}




