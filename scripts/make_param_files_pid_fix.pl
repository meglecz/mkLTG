use warnings;
use strict;
use Data::Dumper;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkLTG;

# Make param files for a combination of different parameters
# pid is fix within each parameter setting

my $outdir = '/home/meglecz/mkLTG/benchmark/param_files_pid_fix/';
my @pid = (100, 97, 95, 90, 85, 80);
my @pcov = (100,90,80,70);
my @phit = (100,90,80,70);

my @taxn = (1,2,3,4,5);
my @refres = (8,7,6,5,4);

my @ltgres = ();
my @seqn = ();

$outdir = add_slash_to_dir($outdir);
my $t = time;
unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}
my $log = $outdir.'make_parameter_space_2.log';
open(LOG, '>', $log);


print LOG "pid\n", Dumper(\@pid), "\n";
print LOG "pcov\n", Dumper(\@pcov), "\n";
print LOG "phit\n", Dumper(\@phit), "\n";
print LOG "taxn\n", Dumper(\@taxn), "\n";
print LOG "seqn\n", Dumper(\@seqn), "\n";
print LOG "refres\n", Dumper(\@refres), "\n";
print LOG "ltgres\n", Dumper(\@ltgres), "\n\n";


my %ind_taxrank = (8, 'species',7,'genus',6,'family',5,'order',4,'class',3,'phylum',2,'kingdom',1,'superkingdom');
my %taxrank_ind = ('species',8,'genus',7,'family',6,'order',5,'class',4,'phylum',3,'kingdom',2,'superkingdom',1);

my $c = 0;
for(my $id = 0; $id < scalar @pid; ++$id)
{
	for(my $cov = 0; $cov < scalar @pcov; ++$cov)
	{
		for(my $hit = 0; $hit < scalar @phit; ++$hit)
		{
			for(my $tn = 0; $tn < scalar @taxn; ++$tn)
			{
				for(my $rr = 0; $rr < scalar @refres; ++$rr)
				{
					++$c;
					my $file = $outdir.'ltg_params_'.$c.'.tsv';
					open(OUT, '>', $file) or die "Cannot open $file\n";
					print OUT "pid	pcov	phit	taxn	seqn	refres	ltgres\n";
					print OUT "$pid[$id]	$pcov[$cov]	$phit[$hit]	$taxn[$tn]	$taxn[$tn]	$ind_taxrank{$refres[$rr]}	species\n";
					close OUT;
				}
			}
		}
	}
}

print LOG "all combinations: $c\n";

print LOG "Runtime ", time - $t, " seconds\n";
close LOG;
exit;

