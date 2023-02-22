use warnings;
use strict;
use Data::Dumper;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkLTG;

# input:
# 	sequences.tsv of COInr
#	taxonomy file of COInr

# ALGO
# 	random select $seqn sequences from COInr (sequences can be restricetd to a given resolution; e.g. get only sequences assigned to species)
#	make a fasta with selected sequences => these will be used later as a query to test mkLTG
#	make a sequence with all COInr sequneces except for the selected ones => reduced_DB
#	make BLASTDB from the reduced_DB
#	BLAST the selected sequences against the reduced database


my %params = 
(
	'db_tsv' => '/home/meglecz/mkCOInr/COInr/COInr.tsv', # seqID	taxID	sequence
	'taxonomy' => '/home/meglecz/mkCOInr/COInr/taxonomy.tsv', # tax_id	parent_tax_id	rank	name_txt	old_tax_id	taxlevel	synonyms
	'seqn' => 1000,
	'select_taxlevel' => 'species', # selet sequences assiged exactly to this resolution
	'outdir' => '/home/meglecz/mkLTG/benchmark/series1',
	# BLAST parameters
	'blast_path' => '',
	'task' => 'megablast',
	'blast_e' => 1e-20,
	'dust' => 'yes',
	'max_target_seqs' => 0, # if zero' =>> 500
	'num_threads' => 8,
	'pid' => 70,
	'qcov_hsp_perc' => 70,
	'delete_tmp' => 1
);
modify_params_from_tags(\%params, \@ARGV);

# cannot be modified from command line
my $outfmt = '6 qseqid sseqid pident length qcovhsp staxids evalue';

my $db_tsv = $params{db_tsv}; 
my $taxonomy = $params{taxonomy}; 
my $seqn = $params{seqn}; 
my $select_taxlevel = $params{select_taxlevel}; 
my $outdir = $params{outdir}; 
# BLAST parameters
my $blast_path = $params{blast_path}; 
my $task = $params{task}; 
my $blast_e = $params{blast_e}; 
my $dust = $params{dust}; 
my $max_target_seqs = $params{max_target_seqs}; 
my $num_threads = $params{num_threads}; 
my $pid = $params{pid}; 
my $qcov_hsp_perc = $params{qcov_hsp_perc}; 
my $delete_tmp = $params{delete_tmp}; 

my %ind_taxrank = (8, 'species',7,'genus',6,'family',5,'order',4,'class',3,'phylum',2,'kingdom',1,'superkingdom');
my %taxrank_ind = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1);


$outdir = add_slash_to_dir($outdir);
unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}
my $date = get_date();
my $t = time;
my $log = $outdir.'make_test_files.log';
open(LOG, '>', $log) or die "Cannot open $log\n";
my @params = print_params_hash_to_log(\%params);
print LOG @params;


####
# read taxonomy to hash
print "Reading taxonomy file\n";
####
my %tax; #$tax{taxid} = (parent_taxid	taxlevel	taxname	taxlevel_index)
my %merged; #$merged{old_taxid} = new_taxid
read_taxonomy_to_hashes($taxonomy, \%tax, \%merged);


####
# random select sequences and make database 
####
	my $fas = $outdir.'random_seq.fas';
	my $db_lo = $outdir.'reduced_DB.fas';
	my $taxids_lo = $outdir.'taxids_reduced_DB.tsv'; # seqID	taxID

if(1)
{
my $select_rank = $taxrank_ind{$select_taxlevel};

	open(IN, $db_tsv) or die "Cannot open $db_tsv\n";
	my @db = <IN>;
	close IN;
	my $title = shift @db; # eliminate title line
	my $n = scalar @db; # total number of sequences


	####
	# random select sequences
	####
	my %random; # $random{random numbers} = ''
	my $c = scalar keys %random;
	while ($c < $seqn)
	{
		my $rn = random_number_1($n);
		unless(exists $random{$rn}) # new random number
		{
			my @tmp = split("\t", $db[$rn]);
			if($tax{$tmp[1]}[3] == $select_rank) # sequence is at a correct taxonomic level
			{
				$random{$rn} = '';
			}
		}
		$c = scalar keys %random;
	}


	####
	# Make fasta file with selected sequences, fasta and taxonomy for making BLAST database with the remaining sequences
	####
	open(FAS, '>', $fas) or die "Cannot open $fas\n";
	open(DBFAS, '>', $db_lo) or die "Cannot open $db_lo\n";
	open(TID, '>', $taxids_lo) or die "Cannot open $taxids_lo\n";
	for(my $i = 0; $i < scalar@db; ++$i)
	{
		my @line = split("\t", $db[$i]);
		if(exists $random{$i})
		{
			print FAS ">$line[0]", '_taxid', $line[1], "\n";
			print FAS $line[2];
		}
		else
		{
			print DBFAS ">$line[0]\n$line[2]";
			print TID "$line[0]\t$line[1]\n";
		}
	}
	close FAS;
	close DBFAS;
	close TID;


	####
	# Make BLAST db and run ltg
	####
	my $makeblastdb = $blast_path.'makeblastdb -dbtype nucl -in '.$db_lo.' -parse_seqids -taxid_map '.$taxids_lo;
	system $makeblastdb;
}
	####
	# BLAST random sequences against LO DB
	####
	my $blastout = $outdir.'blastout_random_seq.tsv';
	local_blast($blast_path, $db_lo, $fas, $blastout, $blast_e, $task, $outfmt, $dust, $qcov_hsp_perc, $pid, $num_threads, $max_target_seqs, 'linux');

if($delete_tmp)
{
	system 'rm '.$db_lo;
	system 'rm '.$taxids_lo;
}

print "Runtime ", time - $t, " seconds\n"; 


exit;

###############################################################
sub random_number_1
{
#If number begins with 0s, zeros are deleted; random number is smaller than number
my ($number) =@_;
	my $digit = length $number;
	my $ok = 1;
	while ($ok == 1)
	{
		my $random_numb = '';
		for (my $i = 0; $i < $digit; ++$i)
		{
			my $l = int rand 10;
#			print "$l ";
		$random_numb .= $l;
 	#end for $i
		}
		if ($random_numb < ($number))
		{
			$random_numb =~ s/^0+//;
			if ($random_numb eq '')
			{
				$random_numb = 0;
			}
			$ok = 0;
			return $random_numb;	
		}
	}
}

########################################################################

sub print_help
{
print '
usage: perl make_LO_db.pl [-options] -in INPUT_FILE -taxonomy TAXONOMY -blast_db BLASTDB
                  -outdir OUTDIR -ltg_params PARMETER_FILE
 arguments:
 
 TO B€ DONE!!!!!!!!!!!!!!!!
 
 
  -in                     name of the input file containg the sequences to be assigned
                          fasta or tsv format (tab separated file with a column titled sequence)
  -taxonomy               tsv file with the following tab separated columns: 
                          taxid parent_taxid taxlevel taxname merged_taxid taxlevel_index
  -ncbitax_dir            directory of ncbi taxonomy dmp files (https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/)
                          only necessary if no -taxonomy file is provided
                          taxonomy file can be created from these files if all sequences in the blast_db have ncbi taxids
  -blast_db               name of the blast database
  -outdir                 name of the otput directory
  -ltg_params             tsv file wih the following tab separated columns:
                          pidentity pcoverage phits min_tax_n min_seq_n min_ref_resolution max_ltg_resolution
  -delete_tmp             0/1; if 1 delete temporary files after the run
 BLAST parameters
  -blast_path             path to blast executables
  -task                   megablast/blastn
  -blast_e                maximum e-value
  -dust                   yes/no
  -max_target_seqs        maximum number of blast hits per query
  -num_threads            number of threads
  -qcov_hsp_perc          minium query coverage in blast
  -pid                    minimum % of identity in blast
  -batch_size             batch size for blast', "\n";
  exit;

}
