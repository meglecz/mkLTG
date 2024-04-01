use warnings;

#################################################
sub add_slash_to_dir
{
 my ($dir, $windows) = @_;

 if($windows)
 {
	 $dir =~ s/\//\\/g;
	 unless($dir eq '')
	 {
		 unless($dir =~  /\\$/)
		 {
			$dir .= '\\';
		 }
	 }
 }
 else
 {
	 unless($dir eq '')
	 {
		 unless($dir =~  /\/$/)
		 {
			$dir .= '/';
		 }
	 }
 }
 return $dir;
}

###############################################################

sub check_and_read_fasta
{
	my ($in) = @_;
	open (IN, $in) or die "Cannot open $in\n";
	
	my $id = '';
	my %seq;
	# reads sequences to hash
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		if($line =~ /^>([^\s]+)/)
		{
			$id= $1;
			if(exists $seq{$id})
			{
				print "ERROR: $id is present more than once in the fasta file\n";
				exit;
			}
		}
		else
		{
			$seq{$id} .= $line;
		}
	}
	
	# check sequence format
	if(scalar keys %seq)
	{
		if(scalar keys %seq == 1 and exists $seq{''}) # wrong format, all lines are read under an id ''
		{
			print "ERROR: input file is either empty or not a correct fasta file\n";
			exit;
		}
		foreach my $id (keys %seq)
		{
			$seq{$id} = uc $seq{$id};
			if($seq{$id} =~ /[^ATCGRYMKBHDVN]/)
			{
				print "ERROR: sequence $id contains unrecognized characters\n";
				exit;
			}
		}
	}
	else
	{
		print "ERROR: input file is either empty or not a correct fasta file\n";
		exit;
	}
	close IN;
	return %seq;

}

###############################################################

sub check_and_read_tsv
{
	my ($in) = @_;
	
	# get the index of the column with the sequences
	open (IN, $in) or die "Cannot open $in\n";
	my $title = <IN>;
	
	$title =~ s/\s*$//;
	$title =~ s/"//g;
	my @title = split("\t", $title);
	my $seq_i = -1; # get the index of the column with the sequences
	for(my $i = 0; $i < scalar @title; ++$i)
	{
		if ($title[$i] =~ /^sequence$/i)
		{
			$seq_i = $i;
			last;
		}
	}
	if($seq_i < 0)
	{
		print "ERROR: Input tsv file does not have a column with 'sequence' as heading\n";
		exit;
	}
	
	
	# read sequences to a hash
	my %seq_inv; # seq_inv{sequence} = id
	while(my $line = <IN>)
	{
		$line =~ s/"//g;
		$line =~ s/\s*$//;
		my @line = split("\t", $line); 
		if(scalar @line <= $seq_i)
		{
			print "ERROR: $line does not have sequence in column ", $seq_i+1, "\n";
			exit;
		}
		else
		{
			$line[$seq_i] = uc $line[$seq_i];
			if($line[$seq_i] =~ /[^ATCGRYMKBHDVN]/) # check sequence format
			{
				print "ERROR: The sequence $line[$seq_i] contains unrecognized characters\n";
				exit;
			}
			else
			{
				++$seq_inv{$line[$seq_i]};
			}
		}
	}
	
	# make %seq using arbitrary numerical seqids
	my %seq;
	my $id = 1;
	if(scalar keys %seq_inv)
	{
		foreach my $s (keys %seq_inv)
		{
			$seq{$id} = $s;
			++$id;
		}
	}
	else
	{
		print "ERROR: There are no correct sequences in the input file in the sequence column\n";
		exit;
	}
	close IN;
	return %seq;

}

###############################################################

sub check_input_file_type_and_read_seqs
{
	my ($in, $seq) = @_;
	
	my $input_format = '';

	if($in =~ /\.fas/i or $in =~ /\.fs/i)
	{
		$input_format = 'fasta'; # get a correct value for $input_format
		%$seq = check_and_read_fasta($in);
	}
	elsif($in =~ /\.[tc]sv/i)
	{
		$input_format = 'tsv';
		%$seq = check_and_read_tsv($in);
	}
	else
	{
		print "The input file should be either a fasta file or a tsv file (tab separated columns and one of the columns must have 'sequence' as a heading)\n";
		exit;
	}
	return $input_format;
}

#################################################

sub delete_file
{
	my ($file, $windows) = @_;
	
#	print "Deleting temporary files\n";
	if($windows)
	{
		system 'del "'.$file.'"';
	}
	else
	{
		system 'rm '.$file;
	}
}


#################################################
sub get_date
{
	my @date = localtime;

	my $y = $date[5] + 1900;
	my $m = $date[4] + 1;
	if((length $m) == 1)
	{
		$m = '0'.$m;
	}
	if((length $date[3]) == 1)
	{
		$date[3] = '0'.$date[3];
	}
	if((length $date[2]) == 1)
	{
		$date[2] = '0'.$date[2];
	}
	if((length $date[1]) == 1)
	{
		$date[1] = '0'.$date[1];
	}
	if((length $date[0]) == 1)
	{
		$date[0] = '0'.$date[0];
	}
	my $d = $y.'-'.$m.'-'.$date[3].'-'.$date[2].'-'.$date[1].'-'.$date[0];
	return $d;
}


#################################################
sub get_file_list_from_folder
{
 my ($folder, $file_motif) = @_;
 
  unless ( opendir(FOLDER, $folder) )
  {
      print "Cannot access to folder $folder\n";
      exit;
  }

my @filenames = grep ( !/^\.\.?$/, readdir(FOLDER) );

#print "@filenames\n";
closedir(FOLDER);
my @files = ();
foreach my $file (sort @filenames)
{
	if ($file =~ /$file_motif/)
	{
		push(@files, $file);
	}
}
@filenames = ();

#print "@files\n";
return @files;
}

###############################################################
sub get_lineage_list
{
	my ($taxid, $tax) = @_;
	#$tax{taxid} = (parent_taxid	taxlevel	taxname	taxlevel_index)

	my @lin; # list of taxid in the linages starting from the largest, ending with the smallest
	while($taxid != 1)
	{
		unshift(@lin, $taxid);
		my $tpar = $$tax{$taxid}[0];
		$taxid = $tpar;
	}
	return @lin;
}

##############################################

sub get_taxlevel
{
	my ($taxid, $tax_ranked_lineage, $tax_rank) = @_;
	#my %tax_ranked_lineage =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
	my $taxlevel = 0;
	
	my %hash = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1); # numerical score for each majot tax level
	my %hash2 = (0=>9, 1=>8, 2 => 7,  3=> 6, 4=> 5, 5=> 4,  6=>3 , 7=> 2, => 8, => 1,  9=> 0); # transform position in @{$ranked_lin{taxid}} => to numerical score

	if(exists $hash{$$tax_rank{$taxid}}) # get numerical score for major taxlevels
	{
		$taxlevel =  $hash{$$tax_rank{$taxid}};
	}
	else # taxlevel between major taxlevels => get an itermediate score (e.g. 7.5 for taxa between genus and species, 8.5 for bellow species)
	{
		for(my $i = 1; $i < 9; ++$i)
		{
			if($$tax_ranked_lineage{$taxid}[$i])
			{
				$taxlevel = $hash2{$i} + 0.5;
				last;
			}
		}
	}
	return $taxlevel;
}

###############################################################
sub local_blast
{
my ($db, $input_fas, $outfile, $e, $task, $outfmt, $dust, $qcov_hsp_perc, $perc_identity, $num_threads, $max_target_seqs) = @_;

	my $blast = 'blastn -task '.$task.' -db "'.$db.'" -query "'.$input_fas.'" -evalue '.$e.' -out "'.$outfile.'" -outfmt "'.$outfmt.'" -dust '.$dust.' -qcov_hsp_perc '.$qcov_hsp_perc.' -perc_identity '.$perc_identity;
	if($num_threads)
	{
		$blast .= ' -num_threads '.$num_threads;
	
	}
	if($max_target_seqs)
	{
		$blast .= ' -max_target_seqs '.$max_target_seqs;
	}
	system $blast;
}

###############################################################

sub ltg
{
	my ($qid, $blastres, $tax, $pid, $pcov, $phit, $taxn, $seqn, $refres, $ltgres) = @_;
	
#	print "$qid, $pid, $pcov, $phit, $taxn, $seqn, $refres, $ltgres\n";

	#$tax{taxid} = (parent_taxid	taxlevel	taxname	taxlevel_index)
	#$blastres{qid} = list of (pident length qcovhsp staxids evalue)
	
	###
	#filter hits by $pid, $pcov, $refres
	###
	# go through all blast hits of the $qid and collect the taxids of the hits that pass the filter
	my %taxids; # $taxids{$staxid} = number of sequneces => get the list of taxids that pass the criteria
	foreach my $line ( @{$$blastres{$qid}}) # $line is a ref to a list, that corresponds to one blast hit cut up to list
	{
#		print $$line[0]; # pident of a blast hit
		if( $$line[0] >= $pid and $$line[2] >= $pcov) # min %identity and %min coverage of the hit is OK
		{
			if( $$tax{$$line[3]}[3] >= $refres) # resolution of the sbj is equal or higher than $refres
			{
				++$taxids{$$line[3]};
			}
		}
	}
	###
	#return 0 if the number of validated sequences (hits) or the number of validated taxa is lower than $seqn or $taxn
	###
	my $hit_n = sum(values %taxids);
	my $tax_n = scalar keys %taxids;
	if($tax_n < $taxn  or $hit_n < $seqn )
	{
		return 0;
	}
	
	###
	# read taxids and lineagaes to a list of hash
	###
	# get the lineages of all taxids, and make a list of hash:
	my @pos_tax; # $pos_tax[taxlevel]{taxid_i} each element is a rank, first is the highest
	# each element is a hash of taxid_i of that level, with the value of number of times the taxid_i is present at that level
	my @taxids = sort {$a <=> $b} keys %taxids;
	for(my $t = 0; $t < $tax_n; ++$t) # go throught all taxids
	{
		my $taxid = $taxids[$t];
		my @lineage = get_lineage_list($taxid, $tax);
		
		for(my $i = 0; $i < scalar @lineage; ++$i)
		{
			$pos_tax[$i]{$lineage[$i]} +=  $taxids{$taxid}; # count the taxid as many time it is present among the hits
		}
	}
#	print Dumper(\@pos_tax);

	
	###
	# get ltg that includes $phit of the hits
	###
	my $ltg_taxid = 0;
	for(my $i = 0; $i < scalar@pos_tax; ++$i)
	{
		my $taxid_i = 0; # local taxid for a given taxlevel
		foreach my $tid ( keys %{$pos_tax[$i]})
		{
			if(($pos_tax[$i]{$tid} / $hit_n) >= $phit) # $tid proportion among the hits is higher than $phit
			{
				$taxid_i = $tid;
				$ltg_taxid = $tid;
				last;
			}
		}
		unless($taxid_i) # no taxid pass the $phit criteria
		{
			last;
		}
	}
	
	###
	# if resolution of the ltg is lower than $ltgres, get the taxon from its lineage that corresponds to $ltgres
	### 
	#$tax{taxid} = (parent_taxid	taxlevel	taxname	taxlevel_index)
	if($ltg_taxid)
	{
		if( $$tax{$ltg_taxid}[3] > $ltgres) # if $ltg_taxid is higher than $ltgres => give the taxid of $ltgres level
		{
			my @lineage = get_lineage_list($ltg_taxid, $tax);
			for(my $i = (scalar @lineage)-2; $i > -1; --$i) # start from the taxid that is a direct parent of the $ltg_taxid
			{
				if($$tax{$lineage[$i]}[3] <= $ltgres)
				{
					$ltg_taxid = $lineage[$i];
					last;
				}
			}
		}
	}
#	print $ltg_taxid, "\n";
	return $ltg_taxid;
	
	

}

###############################################################

sub makedir
{
	my ($dir, $windows) = @_;
	
	
	unless(-e $dir)
	{
		if($windows)
		{
			system 'mkdir "'.$dir.'"';
		}
		else
		{
			system 'mkdir -p "'.$dir.'"';
		}
	}
}


###############################################################

sub make_fastas
{
	my ($seq, $batch_size, $tmpdir) = @_;

	my @ids = sort keys %$seq;
	
	my @fastas;
	for(my $i = 0; $i < scalar @ids; $i = $i + $batch_size)
	{
		my $fas =  $tmpdir.'tempfas_'.$i.'.fasta';
		push(@fastas, $fas);
		open(OUT, '>', $fas) or die "Cannot open $fas\n";
		for(my $j = $i; $j < $i + $batch_size and $j < scalar @ids; ++$j)
		{
			print OUT ">$ids[$j]\n$$seq{$ids[$j]}\n";
		}
		close OUT;
	}
	return @fastas;
}

###############################################################
sub make_ranked_lineage
{
	my ($taxid, $tax, $ind_taxrank) = @_;
#my %ind_taxrank = (8, 'species',7,'genus',6,'family',5,'order',4,'class',3,'phylum',2,'kingdom',1,'superkingdom');

	my @lineage = get_lineage_list($taxid, $tax);
	shift(@lineage); # delete cellular organism
	
	my %ranks;
	my @ranked_lin = ('','','','','','','','');
	foreach my $tid (@lineage) # add taxa in the main levels
	{
		unless($$tax{$tid}[3] =~ /\./) # taxlevel index is an integer
		{
			$ranked_lin[$$tax{$tid}[3] - 1] = $$tax{$tid}[2];
		}
	}

	my $lower_taxa = $ranked_lin[7];
	for(my $i = 6; $i > -1; --$i) # complete lineage if missing high level taxa
	{
		my $taxlevel = $i +1;
		if($ranked_lin[$i] eq '' and $lower_taxa ne '')
		{
			$ranked_lin[$i] = $$ind_taxrank{$taxlevel}.'_'.$lower_taxa;
		}
		$lower_taxa = $ranked_lin[$i];
	}
	return @ranked_lin;
}



##############################################
sub make_taxonomy_with_rank_levels
{
	my ($ncbitax_dir, $outdir, $date) = @_;
	
	# output is a taxonomy file with the following columns:
	# tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel
	# one extra line per merged taxids (old_tax_id)
	
	my $ncbi_tax_rankedlin = $ncbitax_dir.'rankedlineage.dmp';
	my $ncbi_tax_names = $ncbitax_dir.'names.dmp';
	my $ncbi_tax_nodes = $ncbitax_dir.'nodes.dmp';
	my $merged_dump = $ncbitax_dir.'merged.dmp';
	my $taxidlineage = $ncbitax_dir.'taxidlineage.dmp';
	
	my %tax_par; # $tax_par{taxid} = taxid parent
	my %tax_rank; # $tax_rank{taxid} = rank
	my %tax_ranked_lineage; # $ranked_lin{taxid} =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
							# taxname is a scietific name of the taxid, there is allways one but only one scientific name for the taxon
	my %merged; # $merged{taxid} = (list of old taxids)

	# $tax_ranked_lin{taxid} =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
	# taxname is a scietific name of the taxid, there is allways one but only one scientific name for the taxon
	# $tax_ranked_lineage{taxid}[0] scientific name of the taxon
	read_ncbitax_rankedlin($ncbi_tax_rankedlin, \%tax_ranked_lineage);
	# $tax_par{taxid} = taxid parent
	# $tax_rank{taxid} = rank
	read_new_nodes_dmp_to_hash($ncbi_tax_nodes, \%tax_par, \%tax_rank);
	# $name_taxids{name} = (list of taxids) # homonymy
	# $taxid_names{taxid} = (list of all synonymes)
#	read_new_names_dmp_to_hash($ncbi_tax_names, \%name_taxids, \%taxid_names);
	# $taxid_full_lineage_id{taxid} = (list of taxids starting from the most distant)
#	read_new_taxidlineage_dmp_to_hash($taxidlineage, \%taxid_full_lineage_id);


	read_new_merged_dmp_to_hash($merged_dump, \%merged);
#	print Dumper(\%tax_par);

	my $taxonomy = $outdir.'taxonomy_'.$date.'.tsv';
#	tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel

	open(OUT, '>', $taxonomy) or die "Cannot open $taxonomy\n";
	print OUT "tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel\n";
	foreach my $taxid (sort {$a <=> $b} keys %tax_rank)
	{
		my $taxlevel = get_taxlevel($taxid, \%tax_ranked_lineage, \%tax_rank);
		if(exists $merged{$taxid})
		{
			foreach my $old_taxid (@{$merged{$taxid}})
			{
				print OUT "$taxid	$tax_par{$taxid}	$tax_rank{$taxid}	$tax_ranked_lineage{$taxid}[0]	$old_taxid	$taxlevel\n";
			}
		}
		else
		{
			print OUT "$taxid	$tax_par{$taxid}	$tax_rank{$taxid}	$tax_ranked_lineage{$taxid}[0]		$taxlevel\n";
		}
	}
	close OUT;
	return $taxonomy
}

#######################################################
sub modify_params_from_tags
{
	my ($param, $inp) = @_;

	my @bad_tags = ();
	my $version = 0;
	my $help = 0;
	for(my $i = 0; $i<scalar@$inp; $i=$i+2)
	{
		$$inp[$i] =~ s/^-*//;
		if($$inp[$i] =~ /version/i)
		{
			$version = 1;
		}
		elsif($$inp[$i] eq 'h' or $$inp[$i] =~ /help/i)
		{
			$help = 1;
		}
		else
		{
			if(exists $$param{$$inp[$i]})
			{
				$$param{$$inp[$i]} = $$inp[$i+1];
			}
			else
			{
				push(@bad_tags, $$inp[$i]);
			}
		}
	}
	if($version)
	{
		print_version();
		exit;
	}
	if(scalar @bad_tags > 0)
	{
		print_help();
		print "The following tags are not accepted: \n", join(' ', @bad_tags), "\n";
		exit;
	}
	if($help)
	{
		print_help();
		exit;
	}

}

###############################################################

sub print_ltg_fasta_input
{
	my ($ltg, $out, $seq) = @_;
	
	open(OUT, '>', $out) or die "Cannot open $out\n";
	print OUT "seqid	pid	ltg_taxid	ltg_name	ltg_rank	superkingdom	kingdom	phylum	class	order	family	genus	species	sequence\n";
	foreach my $id (sort keys %$seq)
	{
		unless(exists $$ltg{$id})
		{
			$$ltg{$id} = "\t" x 11;
		}
		print OUT $id, "\t", $$ltg{$id}, "\t$$seq{$id}\n";
	}

	close OUT;
}

###############################################################

sub print_ltg_tsv_input
{
	my ($ltg, $out, $seq, $in) = @_;
	# $ltg{seqid} = $pid	$ltg_taxid	$tax{$ltg_taxid}[2]	$tax{$ltg_taxid}[1]	", join("\t", @ranked_lin)
	# $seq{seqid} = seq;
	
	#$seq_inv{sequence} = id
	my %seq_inv = return_hash($seq);
	
	# get the title line and find sequence column index
	open(IN, $in) or die "Cannot open $in\n";
	open(OUT, '>', $out) or die "Cannot open $out\n";
	my $title = <IN>;
	$title =~ s/\s*$//g;
	$title =~ s/"//g;
	my @title = split("\t", $title);
	#and find sequence column index
	my $seq_i = -1; # get the index of the column with the sequences
	for(my $i = 0; $i < scalar @title; ++$i)
	{
		if ($title[$i] =~ /^sequence$/i)
		{
			$seq_i = $i;
			last;
		}
	}
	# complete the title line with ltg headings
	splice(@title, $seq_i, 0, ('pid', 'ltg_taxid', 'ltg_name', 'ltg_rank', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'));
	print OUT join("\t", @title), "\n";
	
	# read $in file
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//g;
		$line =~ s/"//g;
		my @line = split("\t", $line);
		$line[$seq_i] = uc $line[$seq_i]; 
		my $id = $seq_inv{$line[$seq_i]}; # get temporary sequence id used during blast
		my $tmp = "\t" x 11; # make a string with taxonomic info
		if(exists $$ltg{$id})
		{
			$tmp = $$ltg{$id};
		}
		splice(@line, $seq_i, 0, $tmp); # add taxonomic info to the line just before the sequence column
		print OUT join("\t", @line), "\n";
	}
	close IN;
	close OUT;
}

#######################################################
sub print_params_hash_to_log
{
	my ($params) = @_;

	my @print;
	push(@print, "\n####\nPARMATERS:\n\n");
	foreach my $param (sort keys %$params)
	{
		push(@print, "$param: $$params{$param}\n");
	}
	push(@print, "\n####\n\n");
	return @print;
}
#######################################################


sub print_version
{
	print "####################\nmkLTG-0.1.2\n";
	print "April 11, 2023\n####################\n";
}

###############################################################

sub read_blast_results_to_hash
{
	my ($blast, $merged) =  @_;
	
	# $blastres{qid} = list of (pident length qcovhsp staxids evalue)
	my %blastres;
	my %tmp; # $tmp{qid}{sid} = '';
		
	open(IN, $blast) or die "Cannot open $blast\n";
	#qseqid sseqid pident length qcovhsp staxids evalue
	while(my $line = <IN>)
	{
		chomp $line;
		my @line = split("\t", $line);
		my $qid = shift(@line);
		my $sid = shift(@line);
		unless(exists $tmp{$qid}{$sid}) # if multiple hits between sequences take only the first
		{
			if(exists $$merged{$line[3]}) # change merged taxids to valid one
			{
				$line[3] = $$merged{$line[3]};
			}
			push(@{$blastres{$qid}}, \@line);
			$tmp{$qid}{$sid} = '';
		}
	}
	close IN;
	return %blastres;
}

###############################################################

sub read_ltg_params_to_hash
{
	my ($file, $ltg_params, $taxrank_ind) = @_;
#	my %taxrank_ind = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1);
#	my %ltg_params; # $ltg_params{pid}{pcov/phit/taxn/seqn/refres/ltgres} = value

# %hash{pcov/phit/taxn/seqn/refres/ltgres/blast_type} = column number
	my %hash=('pid', '', 'pcov', '', 'phit', '', 'taxn', '', 'seqn', '', 'refres', '', 'ltgres', '');

# check the title line if all parmatares are present and correctly spelled
	open(IN, $file) or die "Cannot open $file\n";
	my $title = <IN>;
	$title =~ s/"//g;
	$title =~ s/\s*$//;
	my @title = split("\t", $title);
	for(my $i = 0; $i < scalar @title; ++$i)
	{
		if(exists $hash{$title[$i]})
		{
			$hash{$title[$i]} = $i; # add column number as value
		}
		else
		{
			print "ERROR: $title[$i] is not recognized as a paramater\n"; # title of the column does not correspond to a valid param
			exit;
		}
	}
	foreach my $param (keys %hash)
	{
		if($hash{$param} eq '') # param is not listed in the input file  
		{
			print "ERROR: $param is not defined in the parameter file $file\n";
			exit;
		}
	}
	
	# read parameter values to hash
	while(my $line = <IN>)
	{
		$line =~ s/"//g;
		$line =~ s/\s*$//;
		my @line = split("\t", $line);
		unless(scalar @title == scalar @line)
		{
			print "ERROR: $file cannot have empty cells\n";
			exit;
		}
		
		foreach my $param (keys %hash)
		{
			my $col_i = $hash{$param};
			my $pid = $hash{pid};
			if($param eq 'ltgres' or $param eq 'refres' ) # change the taxlevels to resulotion index in @line
			{
				$line[$col_i] = $$taxrank_ind{$line[$col_i]};
			}
			elsif($param eq 'phit') # transform %of hits to proportion
			{
				$line[$col_i] = $line[$col_i]/ 100;
			}
			unless($param eq 'pid') # fill out %ltg_params
			{
				$$ltg_params{$line[$pid]}{$param} = $line[$col_i];
			}
		}
	}
	close IN;
	
	my @pident = sort {$a <=> $b} keys %$ltg_params;
	return $pident[0]; # return the lowest value of %identity => this will be used in blast to avoid hits with lower values
}

#########################################################
sub read_ncbitax_rankedlin
{
my ($file, $ranked_lin) = @_;
# $ranked_lin{taxid} =  (tax_name,species,genus,family,order,class,phylum,kingdom,superkingdom)
# taxname is a scietific name of the taxid, there is allways one but only one scientific name for the taxon

	unless(open(IN, $file))
	{
		print "Cannot open $file\n";
	}

	while(my $line = <IN>)
	{
		chomp $line;
		$line =~ s/\s$//;
		$line =~  s/\t\|//g;
		my @line = split("\t", $line);
		my $taxid = shift@line;
		@{$$ranked_lin{$taxid}} = @line;
	}
	close IN;
}

#############################################
sub read_new_nodes_dmp_to_hash
{
my ($file, $par, $rank) = @_;

	# $tax_par{taxid} = taxid parent
	# $tax_par{taxid} = rank
	
unless(open(IN, $file))
{
	print "Cannot open $file\n";
}

while(my $line = <IN>)
{
	chomp $line;
	$line =~ s/\s$//;
	$line =~  s/\t\|//g;
	my @line = split('\t', $line);
	
	$$par{$line[0]} = $line[1];
	$$rank{$line[0]} = $line[2];
}
close IN;
}

#############################################
sub read_new_merged_dmp_to_hash
{
my ($file, $merged) = @_;
	# $merged{valid_taxid} = (old taxid list)


	open(IN, $file) or die "Cannot open $file\n";


	while(my $line = <IN>)
	{
		chomp $line;
		$line =~ s/\s$//;
		$line =~  s/\t\|//g;
		my @line = split('\t', $line);;
		push(@{$$merged{$line[1]}}, $line[0]);
	}
	close IN;
}
###############################################################
sub read_taxonomy_to_hashes
{
	my ($taxfile, $tax, $merged) = @_;
	# taxfile: 
	#taxid	parent_taxid	taxlevel	taxname	merged_taxid	taxlevel_index
	#-1000	-997	species	Poachelas striatus		8
	# $tax{taxid} = (parent_taxid	taxlevel	taxname	taxlevel_index)
	# $merged{old_taxid} = new_taxid
	
	open(IN, $taxfile) or die "Cannot open $taxfile\n";
	my $title = <IN>;
	while(my $line = <IN>)
	{
		chomp $line;
		my @line = split("\t", $line);
		my $taxid = shift(@line);
		my $merged_taxid = splice(@line, 3, 1);
		@{$$tax{$taxid}} = @line;
		if($merged_taxid)
		{
			$$merged{$merged_taxid} = $taxid;
		}
	}
	close IN;
}

###############################################################

sub return_hash
{
	my ($seq) = @_;
	
	my %hash;
	
	foreach my $id (keys %$seq)
	{
		$hash{$$seq{$id}} = $id;
	}
	return %hash;
}

###############################################################
sub sum
{
	my $s = 0;
	foreach my $e (@_)
	{
		$s += $e;
	}
	return $s;
}


1;
