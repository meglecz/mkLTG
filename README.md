# mkLTG
**Taxonomic assignment of metabarcoding sequences using variable % identity thresholds**

## Install

**Perl interpreter** 
This is generally installed in all unix like systems. Otherwise see https://www.activestate.com/products/perl/

**ncbi-blast+**
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/


**mkLTG**
Clone the github repository
~~~~
git clone https://github.com/meglecz/mkLTG.git
~~~~

If you are not comfortable with git, just download the zipped package (**Code** at the top right of the github page https://github.com/meglecz/mkLTG, then **Download ZIP**), dezip the file and you are ready.



## Reference database and taxonomy files
A BLAST database can be formated from a fasta file and a taxmap file
~~~~
makeblastdb -dbtype nucl -in FASTA_FILE -parse_seqids -taxid_map TAXMAP_FILE
~~~~

Format of the **taxmap file** (Sequence identifiers and taxIDs separated by tabulation):

~~~~
HQ449147_1	1000596
6930156	1000596
HQ449148_1	1000597
...
~~~~



**COI database**

The COInr database is a comprehensive non-redundant database with COI sequences from NCBI-nt and BOLD databases (https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13756).

COInr already formatted for mkLTG (and VTAM) is available here:  https://osf.io/vrfwz/

Alternatively, you can download the COInr database from zenodo https://zenodo.org/record/6555985, and prepare a custom database using mkCOInr (https://github.com/meglecz/mkCOInr)



## mkLTG

### Test files
For the sake of this tutorial I have the following file system:

~~~~
~~~~

The files in COInr_for_vtam_2022_05_06_dbV5 were downloaded form https://osf.io/vrfwz/

All other files are provided with mkLTG.

### Run
**Unix like systems**

~~~
cd mkLTG

perl scripts/mkLTG.pl -in data/test1.fas -taxonomy COInr_for_vtam_2022_05_06_dbV5/COInr_for_vtam_taxonomy.tsv -blast_db COInr_for_vtam_2022_05_06_dbV5/COInr_for_vtam -outdir out -out_name test -ltg_params params/params_default.tsv
~~~

**Windows**

~~~
cd mkLTG

perl scripts\mkLTG.pl -in data\test1.fas -taxonomy COInr_for_vtam_2022_05_06_dbV5\COInr_for_vtam_taxonomy.tsv -ncbitax_dir '' -blast_db COInr_for_vtam_2022_05_06_dbV5\COInr_for_vtam -outdir out -out_name test -ltg_params params\params_default.tsv -windows 1
~~~

### Algorithm

Table 1 Default parameter setting

| pid  | pcov | phit | taxn | seqn | refres  | ltgres  |
| ---- | ---- | ---- | ---- | ---- | ------- | ------- |
| 100  | 70   | 70   | 1    | 1    | species | species |
| 97   | 70   | 70   | 1    | 1    | species | species |
| 95   | 70   | 70   | 2    | 2    | species | species |
| 90   | 70   | 70   | 3    | 3    | genus   | species |
| 85   | 70   | 70   | 4    | 4    | family  | genus   |
| 80   | 70   | 70   | 4    | 4    | family  | genus   |

Query sequences are BLASTed against the reference database. Hits are then validated if 
 - % of identity is >= **pid**
 - %of coverage is >=**pcov**
 - the resolution of the reference sequence is at least **refres**

The Lowest Taxonomic Group (**LTG**) is determined if among the validated hits there are at least **taxn** different taxa and **segn** sequences. LTG must contain at least **phit** proportion of the validated hits.

mkLTG starts at the highest pid value of the parameter setting. If LTG cannot be established using this pid and the associated values of the other parameters (e.g. not enough valid hits) the procedure is reiterated using decreasing pid values till LTG is established, or no more pid is provided in the parameter file.

### Arguments and options

 - **in**: Name of the input file containing the sequences to be assigned. Can be fasta or tsv format (tab separated file with a column titled sequence); Compulsory
 - **blast_db**: name of the blast database;  Compulsory
 - **taxonomy**: TSV file with the following tab separated columns: taxid, parent_taxid, taxlevel, taxname, merged_taxid, taxlevel_index;  Either taxonomy or ncbitax_dir is compulsory
 - **ncbitax_dir**: Directory of ncbi taxonomy dmp files (downloaded from https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/) Only necessary if no **taxonomy** file is provided. The taxonomy file is automatically created from the dmp files and cen be used if all sequences have NCBI taxIDs;  Either taxonomy or ncbitax_dir is compulsory
 - **ltg_params**: tsv file wih the following tab separated columns: pid pcov phit taxn seqn refres ltgres ;  Compulsory
 - **outdir**: name of the otput directory;  Compulsory
 - **out_name**: alpha-numeric string for naming output files;  Optional
 - **blastout**: If empty, run blast; otherwise this file is used to make LTG ;  Optional
 - **delete_tmp**: 0/1; if 1 delete temporary files after the run ;  Optional, 1 by default

BLAST parameters
 - **task**: megablast/blastn ;  Optional, megablast by default
 - **blast_e**: E-value threshold;  Optional, 1e-20 by default
 - **dust**: yes/no; apply dust filters for blast;  Optional, yes by default
 - **max_target_seqs**: maximum number of blast hits per query ;  Optional, 500 by default
 - **num_threads**: number of threads when running blast ;  Optional, 8 by default
 - **qcov_hsp_perc**: minium query coverage in blast ;  Optional, 70 by default
 - **batch_size**: batch size for blast (to avoid too large output files) ;  Optional, 1000 by default
 - **windows**: 1/0 set to one if running on windows;  Optional, 0 by default



## Optimize parameters

### Make Test_Files

Randomly select 1000 sequences (seqn) from the COInr.tsv (db_tsv) file (downloaded from https://zenodo.org/record/6555985). Select only sequences assigned to species (select_taxlevel).

The randomy selected sequences (random_seq.fas) will be used as queries for the evaluation of the parameter settings, and the remaining sequences are formatted for blast (reduced_DB.fas). The random_seq.fas is BLASTed against the reduced_DB.fas database.

~~~~
perl scripts\make_test_files.pl -db_tsv COInr\COInr_10000.tsv -taxonomy COInr\taxonomy.tsv -seqn 100 -select_taxlevel species -outdir optimize\series1 -windows 1 
~~~~

### Make parameter files using variable pid

Make parameter files for a large number of parameter value combinations. The possibe parameter values can be set by editing lines 16-22  in the script. This script produces parameter files with multiple pid values. To reduce the number of parameter files, pcov and phits are fixed within each parameter settings (same values for different pids). seqn is the same as taxn, and ltgres is a resolution higher than refres.

~~~~
perl scripts\make_param_files_pid_variable.pl -windows 1 -outdir optimize\param_files_pid_variable
~~~~

###  Make parameter files using a singe pid per parametre setting

Make parameter files for a large number of parameter value combinations. The possibe parameter values can be set by editing lines 15-19  in the script. This script produces parameter files with a single pid values. To reduce the number of parameter files,  seqn is the same as taxn, and ltgres is a resolution higher than refres.

~~~~
perl scripts\make_param_files_pid_fix.pl -windows 1 -outdir optimize\param_files_pid_fix
~~~~

### Run mkLTG for all parameter files 

Run mkLTG for each of the above produced parameter settings using random_seq.fas as a query and the results of the BLAST against reduced_DB.fas. 

At each major taxonomic level, count the number of sequences correctly assigned (TP), incorrectly assigned (FP) and non-assigned (FN)

#### Variable pid 
~~~~
perl scripts\run_ltg_with_multiple_param_setting.pl -fas optimize\series1\random_seq.fas -taxonomy COInr\taxonomy.tsv -blastout optimize\series1\blastout_random_seq.tsv -outdir optimize\series1\LTG_pid_variable -ltg_params_dir optimize\param_files_pid_variable -windows 1
~~~~

#### Fix pid
~~~~
perl scripts\run_ltg_with_multiple_param_setting.pl -fas optimize\series1\random_seq.fas -taxonomy COInr\taxonomy.tsv -blastout optimize\series1\blastout_random_seq.tsv -outdir optimize\series1\LTG_pid_fix -ltg_params_dir optimize\param_files_pid_fix -windows 1
~~~~

### Calculate F1 score, get the parameter setting with the highest F1, make graphs

Adapt the **compair_params.R** script to your use

