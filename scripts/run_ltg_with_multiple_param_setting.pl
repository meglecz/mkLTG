use warnings;
use strict;
use Data::Dumper;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkLTG;

# input;
# 	different LTG param files
#	output of blast


# ALGO
#	run ltg on the same fasta file (containing known sequences) using multiple parameter settings
# 	compare the ltg with the taxa defined by the taxid of the sequence 

my %params = 
(
	'fas' => '',
	'taxonomy' => '', # tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel	synonyms
	'blastout' => '',
	'outdir' => '',
	'ltg_params_dir' => '',
	'delete_tmp' => 1,
	'windows' => 0
);

modify_params_from_tags(\%params, \@ARGV);

my $fas = $params{fas}; 
my $taxonomy = $params{taxonomy}; 
my $blastout = $params{blastout}; 
my $outdir = $params{outdir}; 
my $ltg_params_dir = $params{ltg_params_dir};
my $delete_tmp = $params{delete_tmp};
my $windows = $params{windows};

my %ind_taxrank = (8, 'species',7,'genus',6,'family',5,'order',4,'class',3,'phylum',2,'kingdom',1,'superkingdom');
my %taxrank_ind = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1);
my @param_list = ('pid','pcov','phit','taxn','seqn','refres','ltgres');

$outdir = add_slash_to_dir($outdir, $windows);
makedir($outdir, $windows);
$ltg_params_dir = add_slash_to_dir($ltg_params_dir, $windows);

my $log = $outdir.'run_ltg_params.log';
open(LOG, '>', $log) or die "Cannot open $log\n";
my @params = print_params_hash_to_log(\%params);
print LOG @params;
my $t = time;
my @ltg_params_files = get_file_list_from_folder($ltg_params_dir, '^ltg_params');

my %seq = read_fasta_to_hash_id_till_first_space($fas);
####
# read taxonomy to hash
print "Reading taxonomy file\n";
####
my %tax; #$tax{taxid} = (parent_taxid	taxlevel	taxname	taxlevel_index)
my %merged; #$merged{old_taxid} = new_taxid
read_taxonomy_to_hashes($taxonomy, \%tax, \%merged);

#print $tax{218371}[3], "\n";

#my $stat_resolution = $outdir.'stat_resolution.tsv';
#open(RES, '>', $stat_resolution) or die "Cannot open $stat_resolution\n";
#print RES "param_setting	taxlevel	resolution_ltg	highest matching level	number of sequences\n";
my $stat_match = $outdir.'stat_match.tsv';
open(MAT, '>', $stat_match) or die "Cannot open $stat_match\n";
print MAT "param_setting	", join("\t", @param_list);
#print MAT "	taxlevel	taxlevel_ind	TP (correct)	FP (incorrect)	FN (missing)	FDR (FP/(FP+TP))	TPR (TP/(TP+FN))	PPV (TP/(TP+FP))	F1 (2*TPR*PPV/(TPR+PPV))	(FP+FN)/All\n";
print MAT "	taxlevel	taxlevel_ind	TP (correct)	FP (incorrect)	FN (missing)\n";


foreach my $params (sort @ltg_params_files)
{
	my $ltg_params = $ltg_params_dir.$params;
	my $out_name = $params;
	$out_name =~ s/ltg_params_//;
	$out_name =~ s/\.tsv//;
	print "\n####\nRunning params $out_name\n####\n";
	
##########################

	my %ltg_params; # $ltg_params{pid}{pcov/phit/taxn/seqn/refres/ltgres} = value
	my %params_print; # $params_print{pcov/phit/taxn/seqn/refres/ltgres} = (list of values for each pid)
	my $perc_identity = read_ltg_params_to_hash_bis($ltg_params, \%ltg_params, \%taxrank_ind, \%params_print);
	# $perc_identity is the lowest value in the $ltg_params file

	my @pident = reverse sort {$a <=> $b} keys %ltg_params;
	my %ltg; # $ltg{seqid} = $pid	$ltg_taxid	$tax{$ltg_taxid}[2]	$tax{$ltg_taxid}[1]	", join("\t", @ranked_lin),

#	print Dumper(\%params_print);

	####
	# read blast results to hash
	####
	# $blastres{qid} = list of (pident length qcovhsp staxids evalue)
	# change merged taxids to valid taxids
	my %blastres = read_blast_results_to_hash($blastout, \%merged);
#	print Dumper(\%blastres);

	####
	# ltg
	####
	foreach my $qid (keys %blastres)
	{
		foreach my $pid (@pident) # start for the highest pindent
		{
			my $ltg_taxid = ltg($qid, \%blastres, \%tax, $pid, $ltg_params{$pid}{pcov}, $ltg_params{$pid}{phit}, $ltg_params{$pid}{taxn}, $ltg_params{$pid}{seqn}, $ltg_params{$pid}{refres}, $ltg_params{$pid}{ltgres});
			if($ltg_taxid) # if taxid was deduced
			{
				my @ranked_lin = make_ranked_lineage($ltg_taxid, \%tax, \%ind_taxrank); # get ranked lineage
				my @line = ($pid, $ltg_taxid, $tax{$ltg_taxid}[2], $tax{$ltg_taxid}[1]); # make a list with ltg info
				splice(@line, scalar @line, 0, @ranked_lin); # add ranked lineage to list
				$ltg{$qid} = join("\t", @line); # put list to the hash
				last;
			}
		}
	}
#	print Dumper(\%ltg);

	####
	# print results to an output file
	####
	my $ltg_file = $outdir.$out_name.'_ltg.tsv';
	print_ltg_fasta_input(\%ltg, $ltg_file, \%seq);

###########################

	####
	# compare lineage of the original taxid and ltg of the selected sequences
	####
	
	my $ltg_comp = $outdir.$out_name.'_ltg_completed.tsv';

	my %count; # $count{ranked_lin_level-ltg_ranked_lin_level-ltg_ranked_lin_level-other_seq_same_taxid}
	open(IN, $ltg_file) or die "Cannot open $ltg_file\n";
	open(LTGC, '>', $ltg_comp) or die "Cannot open $ltg_comp\n";

	my $title = <IN>;
	my @title = split("\t", $title);
	my @tl = ('superkingdom_taxid','kingdom_taxid','phylum_taxid','class_taxid','order_taxid','family_taxid','genus_taxid','species_taxid');
	splice(@title, 13, 0, @tl);
	splice(@title, 21, 0, ('resolution_taxid','resolution_ltg','highest matching level'));
	print LTGC join("\t", @title);

	# $comp{taxlevel} = (correct, incorrect, missing)
	my %comp;
	@{$comp{1}} = (0,0,0);
	@{$comp{2}} = (0,0,0);
	@{$comp{3}} = (0,0,0);
	@{$comp{4}} = (0,0,0);
	@{$comp{5}} = (0,0,0);
	@{$comp{6}} = (0,0,0);
	@{$comp{7}} = (0,0,0);
	@{$comp{8}} = (0,0,0);

	while(my $line = <IN>)
	{
		chomp $line;
		my @line = split("\t", $line);
		my @ltg_ranked_lin = splice(@line, 5, 8); # get ranked lineage of ltg
		splice(@line, 5, 0, @ltg_ranked_lin); # put back ranked lineage of ltg
		
		my @id = split('_taxid', $line[0]);
		my $sid = $id[0];
		my $taxid = $id[1];
		my @ranked_lin = make_ranked_lineage($taxid, \%tax, \%ind_taxrank); # get ranked lineage
		
		splice(@line, 13, 0, @ranked_lin); # add anked lineage of taxid

		my $last_match_level = 8; # Highest resolution with matching taxa
		my $ranked_lin_level = 8; # Resolution of taxid based on ranked lineage (ignore other taxlevels)
		my $ltg_ranked_lin_level = 8; # Resolution of ltg  based on ranked lineage (ignore other taxlevels)

		for(my $i = 7; $i >= 0; --$i)
		{
			unless($ranked_lin[$i])
			{
				$ranked_lin_level = $i; # get the resolution of taxid
			}
			if($ltg_ranked_lin[$i]) # taxassign at level $i
			{
				if($ltg_ranked_lin[$i] eq $ranked_lin[$i])
				{
					++$comp{$i+1}[0];
				}
				else
				{
					++$comp{$i+1}[1];
				}
			}
			else # resolution is lower then i
			{
				++$comp{$i+1}[2];
				$ltg_ranked_lin_level = $i; # get the resolution of ltg
			}
			unless($ranked_lin[$i] eq $ltg_ranked_lin[$i])
			{
				$last_match_level = $i;
			}
		}

		my $t = $ranked_lin_level."\t".$ltg_ranked_lin_level."\t".$last_match_level;
		++$count{$t};
		
		splice(@line, 21, 0, ($ranked_lin_level,$ltg_ranked_lin_level,$last_match_level));
		print LTGC join("\t", @line), "\n"; 
	}
	close IN;
	close LTGC;


#	foreach my $c (sort keys %count)
#	{
#		print RES $out_name, "\t", $c, "\t", $count{$c}, "\n";
#	}

	foreach my $tl (sort keys %comp)
	{
		my $fdp = 0; # FP/(FP+TP)
		my $tpr = 0; # TP/(TP+FN)
		my $ppv = 0; # TP/(TP+FP)
		my $f1 = 0; # 2*TPR*PPV/(TPR+PPV)
		my $false = 0;
		if(0)
		{
			if($comp{$tl}[1]+$comp{$tl}[0])
			{
				$fdp = $comp{$tl}[1]/($comp{$tl}[1]+$comp{$tl}[0]); # FP/(FP+TP)
			}
			if($comp{$tl}[0]+$comp{$tl}[2])
			{
				$tpr = $comp{$tl}[0]/($comp{$tl}[0]+$comp{$tl}[2]); # TP/(TP+FN)
			}
			if($comp{$tl}[0]+$comp{$tl}[1])
			{
				$ppv = $comp{$tl}[0]/($comp{$tl}[0]+$comp{$tl}[1]); # TP/(TP+FP)
			}
			if($tpr+$ppv)
			{
				$f1 = 2*$tpr*$ppv/($tpr+$ppv);# 2*TPR*PPV/(TPR+PPV)
			}
			$false = ($comp{$tl}[1]+$comp{$tl}[2])/($comp{$tl}[0]+$comp{$tl}[1]+$comp{$tl}[2]);
		}
		print MAT  $out_name;
		foreach my $par (@param_list)
		{
			print MAT "\t", join(';', @{$params_print{$par}}),
		}
#		print MAT  "\t", $ind_taxrank{$tl}, "\t", $tl, "\t", join("\t", @{$comp{$tl}}),  "\t", $fdp, "\t", $tpr,  "\t", $ppv, "\t", $f1, "\t", $false, "\n";
		print MAT  "\t", $ind_taxrank{$tl}, "\t", $tl, "\t", join("\t", @{$comp{$tl}}),  "\n";

	}

	if($delete_tmp)
	{
		delete_file($ltg_file, $windows);
		delete_file($ltg_comp, $windows);
	}
}# end foreach params
#close RES;
close MAT;

print "Runtime: ", time - $t,  " seconds\n";


exit;



###############################################################

sub read_ltg_params_to_hash_bis
{
	my ($file, $ltg_params, $taxrank_ind, $params_print) = @_;
#	my %taxrank_ind = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1);
#	my %ltg_params; # $ltg_params{pid}{pcov/phit/taxn/seqn/refres/ltgres} = value
# $params_print{pcov/phit/taxn/seqn/refres/ltgres} = (list of values for each pid)

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
				push(@{$$params_print{$param}}, $line[$col_i]);
			}
			elsif($param eq 'phit') # transform %of hits to proportion
			{
				push(@{$$params_print{$param}}, $line[$col_i]);
				$line[$col_i] = $line[$col_i]/ 100;
			}
			else
			{
				push(@{$$params_print{$param}}, $line[$col_i]);
			}
			unless($param eq 'pid') # fill out %ltg_params
			{
				$$ltg_params{$line[$pid]}{$param} = $line[$col_i];
			}
			# fill params_print
			
		}
	}
	close IN;
	
	my @pident = sort {$a <=> $b} keys %$ltg_params;
	return $pident[0]; # return the lowest value of %identity => this will be used in blast to avoid hits with lower values
}

######################################################
sub read_fasta_to_hash_id_till_first_space
{
my ($file) = @_;

open(DATA, $file) or die(print "Cannot open $file\n");

my $code = '';
my %hash = ();
while(my $line= <DATA>)
{	
	chomp $line;
	$line =~ s/\s*$//;
	if ($line =~ />/)
	{
		$code = $line;
		$code =~ s/>//;
		$code =~ s/\s.*//;
	}
	else
	{
		$line =~ s/\s//g;
		$hash{$code} .=$line;
	}
}
close DATA;
 return %hash;
 
}

#######################################################
sub print_help
{

print '
usage: perl ltg.pl [-options] -in INPUT_FILE -taxonomy TAXONOMY -blast_db BLASTDB
                  -outdir OUTDIR -ltg_params PARMETER_FILE
 arguments:
  -fas                    fasta file with query sequences
  -blastout               output of BLAST quary fasta aginst the database without the query sequences
  -taxonomy               tsv file with the following tab separated columns: 
                          taxid parent_taxid taxlevel taxname merged_taxid taxlevel_index
  -outdir                 name of the otput directory
  -ltg_params_dir         directory containing the parameter files: ech of them should have the following columns
                          pid pcov phit taxn seqn refres ltgres
  -delete_tmp             0/1; if 1 delete temporary files after the run
  -windows                set to one if running on windows', "\n";
  exit;


}
