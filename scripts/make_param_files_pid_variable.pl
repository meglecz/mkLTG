use warnings;
use strict;
use Data::Dumper;
use FindBin qw( $RealBin );
use lib "$RealBin";
use mkLTG;

# Make param files for a combination of different parameters
# pid varies in within each parameter setting

 

##########################
##########################
# set the parameter space here
my @pid = (100, 97, 95, 90, 85, 80); # list of pid values; same in all parametre setting
my @pcov = ([100,100,100,100,100,100], [90,90,90,90,90,90], [80,80,80,80,80,80], [70,70,70,70,70,70]); # pcov values for each pid values
my @phit = ([100,100,100,100,100,100], [90,90,90,90,90,90], [80,80,80,80,80,80], [70,70,70,70,70,70]); # phit values for each pid values

#my @taxn = ([1], [1], [2], [3], [4], [4]);
my @taxn = ([1], [1,2], [2,3,4], [3,4], [4,5], [4,5]); # taxn values for each pid values
my @refres = ([8], [8,7], [8,7], [8,7,6], [7,6,5], [6,5,4]); # refres values for each pid values; 8 species, 7 genus, 6 familly, 6 order, 5 class
#my @refres = ([8], [8], [8], [8], [7], [6]);
##########################
##########################
my @ltgres = ();
my @seqn = ();


my %params = 
(
'windows' => 0,
'outdir' => ''
);
modify_params_from_tags(\%params, \@ARGV);

my $windows = $params{windows}; 
my $outdir = $params{outdir};

$outdir = add_slash_to_dir($outdir, $windows);
my $t = time;
makedir($outdir, $windows);

my $log = $outdir.'make_files_pid_variable.log';
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

@taxn = get_all_combinations(\@taxn, 1);
@refres = get_all_combinations(\@refres, -1); 

my $nparam_pcov = scalar @pcov;
my $nparam_phit = scalar @phit;
my $nparam_taxn = scalar @taxn;
my $nparam_refres = scalar @refres;

my $all = $nparam_pcov * $nparam_phit *$nparam_taxn * $nparam_refres;

print LOG "nparam_pcov: $nparam_pcov\n";
print LOG "nparam_phit: $nparam_phit\n";
print LOG "nparam_taxn: $nparam_taxn\n";
print LOG "nparam_refres: $nparam_refres\n";
print LOG "all combinations: $all\n";

my $c = 0;
for(my $pc = 0; $pc < scalar @pcov; ++$pc)
{
	for(my $ph = 0; $ph < scalar @phit; ++$ph)
	{
		for(my $mt = 0; $mt < scalar @taxn; ++$mt)
		{
			for(my $mr = 0; $mr < scalar @refres; ++$mr)
			{
				++$c;
				my $file = $outdir.'ltg_params_'.$c.'.tsv';
				open(OUT, '>', $file) or die "Cannot open $file\n";
				print OUT "pid	pcov	phit	taxn	seqn	refres	ltgres\n";
				for(my $pi = 0; $pi < scalar @pid; ++$pi)
				{
					my $max_ltg_res;
					if($refres[$mr][$pi] == 8)
					{
						$max_ltg_res = 8;
					}
					else
					{
						$max_ltg_res= $refres[$mr][$pi] +1;
					}
					print OUT "$pid[$pi]	$pcov[$pc][$pi]	$phit[$ph][$pi]	$taxn[$mt][$pi]	$taxn[$mt][$pi]	$ind_taxrank{$refres[$mr][$pi]}	$ind_taxrank{$max_ltg_res}\n";
				}
				close OUT
			}
		}
	}
}

print LOG "Runtime ", time - $t, " seconds\n";
close LOG;
exit;


########################################################"

sub get_all_combinations
{
	my ($param_space, $order) = @_;
	# $order -1, keep only parm list with decreasing values, 0 no selection, 1, keep only parm list with increasing values
	
	my %hash;
	# initialize hash
	for(my $j = 0; $j < scalar  @{$$param_space[0]}; ++$j)
	{
		$hash{$$param_space[0][$j]} = '';
	}
	
	# compelte strings
	for(my $i = 1; $i < scalar @$param_space; ++$i)
	{
		foreach my $k (keys %hash)
		{ 
			for(my $j = 0; $j < scalar  @{$$param_space[$i]}; ++$j)
			{
				my $tmp = $k.';'.$$param_space[$i][$j];
				$hash{$tmp} = '';
			}
		}
	}
#	print Dumper(\%hash);

	# delete strigs shorther than scalar @refres
	my @params;
	foreach my $k (sort keys %hash)
	{ 
		my @l = split(';', $k);
		if(scalar @l < scalar @$param_space)
		{
			delete $hash{$k};
		}
		else # list complet
		{
			
			if($order == -1) # decreasing order
			{
				for(my $i = 1; $i < scalar @l; ++$i)
				{
					if($l[$i]>$l[$i-1])
					{
						delete $hash{$k};
						$i = scalar @l;
					}
				}
			}
			if($order == 1) # increasing order
			{
				for(my $i = 1; $i < scalar @l; ++$i)
				{
					if($l[$i]<$l[$i-1])
					{
						delete $hash{$k};
						$i = scalar @l;
					}
				}
			}
			if(exists $hash{$k}) # $k has not been deleted
			{
				push(@params, \@l);
			}
		}
	}
	return @params;
}

########################################################"

sub combine
{
	my ($list_small, $list_big, $restrict) = @_;
	
	my @combined_list; # combined_list{combinations of 2 params} = [[param1], [$param2]]
	for(my $s = 0; $s < scalar @$list_small; ++$s)
	{
		for(my $b = 0; $b < scalar @$list_big; ++$b)
		{
			my $bool =1;
			if($restrict)
			{
				for(my $i = 0; $i < scalar @{$$list_small[$s]}; ++$i)
				{
					if($$list_small[$s][$i] > $$list_big[$b][$i])
					{
						$bool = 0;
						last;
					}
				}
			}
			
			if($bool)
			{
				push(@combined_list, [@{$$list_small[$s]}, @{$$list_big[$b]}]);
			}
		}
	}
	return @combined_list;
}

########################################################################

sub print_help
{
print '
usage: perl make_param_files_pid_variable.pl [-options] -outdir OUTDIR 

  -outdir                 name of the otput directory
  -windows                set to one if running on windows', "\n";
  exit;

}


