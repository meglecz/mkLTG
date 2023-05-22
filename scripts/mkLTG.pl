use warnings;
use strict;
use Data::Dumper;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkLTG;

# conda activate vtam
# Nucleotide BLAST 2.10.1+

# INPUT:
# -fasta file with sequences to be assigned OR tsv file with 'sequence' as heading for the column that contains the sequences to be assigned
# -blast database 
# -taxonomy file : tsv file with the following tab separated columns: taxid parent_taxid taxlevel taxname merged_taxid taxlevel_index
# 	if the taxonomy file does not exists, it can be created based on the ncbi taxonomy dmp files downloaded from
#	https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/. In that case, all sequences in the blast_db must have ncbi taxids
# -ltg_params.tsv : ltg parameters applided for each %identity level
#	tab separated file with the following columns:  pid, pcov, phit, taxn, seqn, refres, ltgres

# ALGO
# make temp fasta file with max $batch_size sequences (to avoid too large blast output files)
# each fasta file is then taxassigned and the results are saved in a hash
# TAXASSIGN
#	start with the highest pid (%identity) => determine ltg if possibe, if not take the next pid
#	ltg => the lowest taxonomic group that contains phit proportion of the validated hits and it should be based on at least taxn taxa and seqn sequences
#		phit: % of hits in the ltg
#		taxn: minimum number of taxa among the validated hits
#		seqn: minimum number of sequences among the validated hits
#	validated hits => pid , pcov , refres are checked
#		pid: minimum % identity between quary and reference sequence
#		pcov: minimum % of query coverage
#		refres: minimun resolution of the reference sequence; 
#			accepted resolutions in increasing order: 
#			root, superkingdom, kingdom, phylum, class, order, family, genus, species

# OUTPUT
# if input was a fasta file => tsv with the following colums:
# seqid	pid	ltg_taxid	ltg_name	ltg_rank	superkingdom	kingdom	phylum	class	order	family	genus	species	sequence
# sequences without ltg are kept in the file

# if input was a tsv file => complete the tsv by adding the following columns before the sequence column:
# pid	ltg_taxid	ltg_name	ltg_rank	superkingdom	kingdom	phylum	class	order	family	genus	species
# sequences without ltg are kept in the file


my %params = 
(
	'in' => '', # fasta or TSV file with sequnces to be assigned
	'taxonomy' => '',# if empty file is created from ncbi tax dmp files
	'ncbitax_dir' => '', # unless $taxonomy, make a taxonomy file including rank levels; can be emty if $taxonomy extists
	'blast_db' => '',
	'outdir' => '',
	'out_name' => '', # alpha-numeric string for naming output files
 	'ltg_params' => '',
	'delete_tmp' => 1,
	'windows' => 0, # set to 1 if running on windows
	'blastout' => '', # if empty, run blast. otherwise this file is used to make ltg
	# BLAST parameters
	'task' => 'megablast',
	'blast_e' => 1e-20,
	'dust' => 'yes',
	'max_target_seqs' => 0, # if zero' =>> 500
	'num_threads' => 8,
	'qcov_hsp_perc' => 70,
	'batch_size' => 1000 # number of sequences in a fasta that are blasted in one go
);

modify_params_from_tags(\%params, \@ARGV);

my $in = $params{in}; 
my $taxonomy = $params{taxonomy}; 
my $ncbitax_dir = $params{ncbitax_dir}; 
my $blast_db = $params{blast_db}; 
my $outdir = $params{outdir}; 
my $out_name = $params{out_name}; 
my $ltg_params = $params{ltg_params};
my $delete_tmp = $params{delete_tmp};  
# BLAST parameters
my $blastout = $params{blastout};  
my $task = $params{task}; 
my $blast_e = $params{blast_e}; 
my $dust = $params{dust}; 
my $max_target_seqs = $params{max_target_seqs}; 
my $num_threads = $params{num_threads}; 
my $qcov_hsp_perc = $params{qcov_hsp_perc}; 
my $batch_size = $params{batch_size}; 
my $windows = $params{windows};

# cannot be modified from command line
my $outfmt = '6 qseqid sseqid pident length qcovhsp staxids evalue';

$outdir = add_slash_to_dir($outdir, $windows);
$ncbitax_dir = add_slash_to_dir($ncbitax_dir, $windows);



my $date = get_date();
my $out = $outdir.'ltg_'.$date.'.tsv';
my $log = $outdir.'ltg_'.$date.'.log';
if($out_name)
{
	$out = $outdir.$out_name.'_ltg.tsv';
	$log = $outdir.$out_name.'_ltg.log';
}

# make temp dir
my $tmpdir = $outdir.'tmp_'.$date;
if($out_name)
{
	$tmpdir = $outdir.$out_name.'_tmp';
}
$tmpdir = add_slash_to_dir($tmpdir, $windows);
makedir($tmpdir, $windows);

open(LOG, '>', $log) or die "Cannot open $log\n";
my @params = print_params_hash_to_log(\%params);
print LOG @params;

my $t = time;

####
# read ltg params to hash
print "Reading ltg parameters\n";
####
my %taxrank_ind = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1);
my %ind_taxrank = (8, 'species',7,'genus',6,'family',5,'order',4,'class',3,'phylum',2,'kingdom',1,'superkingdom');
my %ltg_params; # $ltg_params{pid}{pcov/phit/taxn/seqn/refres/ltgres} = value
# $perc_identity is the lowest value in the $ltg_params file
my $perc_identity = read_ltg_params_to_hash($ltg_params, \%ltg_params, \%taxrank_ind);
print_param_setting_to_log(\%ltg_params);

####
# read input files and check format
print "Reading sequences\n";
####
# the input format is guessed from the extention
my %seq;
my $input_format = check_input_file_type_and_read_seqs($in, \%seq);
# make fasta files by batch of $batch_size
my @fastas = make_fastas(\%seq, $batch_size, $tmpdir);


####
# Make taxonomy file if it does not exists
####
unless($taxonomy)
{
	print "Make taxonomy file from ncbi tax dmp files\n";
	unless(-e $ncbitax_dir)
	{
		print "Taxonomic information should be provided either in -taxonomy or -ncbitax_dir parameters\n";
		exit;
	}
	$taxonomy = make_taxonomy_with_rank_levels($ncbitax_dir, $outdir, $date);
}

####
# read taxonomy to hash
print "Reading taxonomy file\n";
####
my %tax; #$tax{taxid} = (parent_taxid	taxlevel	taxname	taxlevel_index)
my %merged; #$merged{old_taxid} = new_taxid
read_taxonomy_to_hashes($taxonomy, \%tax, \%merged);

####
# make ltg for all sequences in all of the temporary fasta files
# keep info in a %hash
print "Making LTG\n";
####
my @pident = reverse sort {$a <=> $b} keys %ltg_params;
my %ltg; # $ltg{seqid} = $pid	$ltg_taxid	$tax{$ltg_taxid}[2]	$tax{$ltg_taxid}[1]	", join("\t", @ranked_lin),
foreach my $fas (@fastas)
{
	####
	# blast
	####
	
	my $blastout_tmp = $blastout;
	unless($blastout) # run blast only if BLAST output is not given in the input
	{
		$blastout_tmp = $fas;
		$blastout_tmp =~ s/\.fasta/_blastout.tsv/;
		local_blast($blast_db, $fas, $blastout_tmp, $blast_e, $task, $outfmt, $dust, $qcov_hsp_perc, $perc_identity, $num_threads, $max_target_seqs);
	}
	####
	# read blast results to hash
	####
	# $blastres{qid} = list of (pident length qcovhsp staxids evalue)
	# change merged taxids to valid taxids
	my %blastres = read_blast_results_to_hash($blastout_tmp, \%merged);


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
}

####
# print results to an output file
# if input was fasta, make a tsv outfile
# if input was a tsv, copy the content in new tsv outfile and complete the tsv with the taxinfo. Add tax info before the sequence column
print "Writing output\n";
####
if($input_format eq 'fasta')
{
	print_ltg_fasta_input(\%ltg, $out, \%seq);
}
else
{
	print_ltg_tsv_input(\%ltg, $out, \%seq, $in);
}


####
# Deleting temporary files
####

if($delete_tmp)
{
	delete_dir($tmpdir, $windows);
}


print "Runtime: ", time - $t, "s\n";
print LOG "Runtime: ", time - $t, "s\n";
close LOG;
exit;

################################################

sub print_param_setting_to_log
{
	my($ltg_params) = @_;
# $ltg_params{pid}{pcov/phit/taxn/seqn/refres/ltgres} = value

	my @p = ('pcov','phit','taxn','seqn','refres','ltgres');
	print LOG "pid	", join("\t", @p), "\n";
	foreach my $pid (sort {$a <=> $b} keys %$ltg_params)
	{
		print LOG "$pid";
		foreach my $p (@p)
		{
			print LOG "	$$ltg_params{$pid}{$p}";
		}
		print LOG "\n",
	}

}
################################################

sub delete_dir
{
	my ($dir, $windows) = @_;
	
	print "Deleting temporary files\n";
	if($windows)
	{
		system 'del "'.$dir.'"*.tsv';
		system 'del "'.$dir.'"*.fasta';
		system 'rmdir "'.$dir.'" \q';
	}
	else
	{
		system 'rm -r "'.$dir.'"';
	}
}

#######################################################
sub print_help
{

print '
usage: perl mkLTG.pl [-options] -in INPUT_FILE -taxonomy TAXONOMY -blast_db BLASTDB
                  -outdir OUTDIR -ltg_params PARMETER_FILE
 arguments:
  -in                     name of the input fasta or TSV file containing the sequences to be assigned; the TSV file must have "sequence" as a heading for one of the columns
  -taxonomy               tsv file with the following tab separated columns: 
                          taxid parent_taxid taxlevel taxname merged_taxid taxlevel_index
  -ncbitax_dir            directory of ncbi taxonomy dmp files (https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/)
                          only necessary if no -taxonomy file is provided
                          taxonomy file can be created from these files if all sequences in the blast_db have ncbi taxids
  -blast_db               name of the blast database
  -ltg_params             tsv file wih the following tab separated columns:
                          pid pcov phit taxn seqn refres ltgres
  -outdir                 name of the otput directory
  -out_name               alpha-numeric string for naming output files
  -blastout               if empty, run blast; otherwise this file is used to make ltg
  -delete_tmp             0/1; if 1 delete temporary files after the run
  -windows                0/1, set to one if running on windows
 BLAST parameters
  -task                   megablast/blastn
  -blast_e                maximum e-value
  -dust                   yes/no
  -max_target_seqs        maximum number of blast hits per query
  -num_threads            number of threads
  -qcov_hsp_perc          minium query coverage in blast
  -batch_size             batch size for blast' , "\n";
}
