#

# perl v5.20.0
# bedtools 2.24.0

######################
## DESCRIPTION
###################### 
## AUTHOR naveed.ishaque@charite.de
## DATE Jan 2019
######################
## This script parses output from multiple SOPHIA structural variation "minEventScore3" files, and create structures that can be plotted circlize
## Technically it should work with any set of SV calls in BED-PE format (where the first 3 cols are the start of the SV, and the next 3 columns and the end of the SV)
## This script expects a folder sctructure created by the softlink project script created by @ishaque
## The PID list and or file list can be forced using --pid_list and --file_list options
######################

######################
## CHANGELOG
######################
## v0.1 24/01/2018 - creation
## v0.2 25/01/2018 - completed
## v0.3 31/01/2018 - un-completed! Synchronising parsing options with the oncoprint script.... perhaps these scripts should be modularised into a workflow
######################

use strict;
use Getopt::Long;
use Cwd qw(cwd);

my @arguments = @ARGV;

my $usage = "

This script creates a table for producing cohort level SV plots

	perl $0
		
		--results_dir     -r   [path: file path to the results_per_pid folder. REQUIRED!]
		--output_prefix   -o   [string: output prefix. REQUIRED!]
		--sv_close        -s   [integer: distance to define close SV. Default 10000]
		--deldup_window   -d   [integer: size of del/dup windows for aggregation. Default 1000000]
		--pids            -p   [file or list: file containing a list of PIDs, or a space separated list of pids]
		--file_list       -f   [list of files to use... over rides PIDs]
		--help            -h   [print usage information]

Example call:

	perl $0 -r /path/tp/hipo_050/whole_genome_sequencing/results_per_pid -o test_H050 -p \"H050-ABCD H050-EFGH H050-IJKL\"

"; 

my $version = "v0.3"; 

######################
## ENV DEFINITIONS
######################

my $bedtools_version = `bedtools -version`;
chomp ($bedtools_version);

die "ERROR :: bedtools version :: \"expected  > v2.24.0\", but found \"$bedtools_version\"\n" unless ($bedtools_version >= "v2.24.0");

warn "\n# ENV :: perl $^V (expected > v5.18.1)\n";
warn "# ENV :: $bedtools_version (expected > v2.24.0)\n";

######################
## FOLDER FILE VARIABLE DEFINITIONS
######################

# These file definitions should not change as sophia only works with hg19 builds!

my $results_dir   = "/icgc/dkfzlsdf/analysis/hipo/my_hipo_project/whole_genome_sequencing/results_per_pid/";
my $chrom_arms    = "/applications/otp/ngs_share_complete-copy/assemblies/hg19_GRCh37_1000genomes/stats/hg19armsizesXY.plain.bed";
my $genome_fai    = "/applications/otp/ngs_share_complete-copy/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef_Phix/hs37d5_PhiX.fa.fai";
my $genes_bed     = "/applications/otp/ngs_share_complete-copy/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_proteinCodingGenes_plain.bed.gz";

# SV file expected folder, preffix and extension

my $pre_sv    = "svs_";
my $ext_sv  = "_filtered_somatic_minEventScore3.tsv"; # comment: default from OTP
my $folder_sv    = "SOPHIA_";

######################
## DEFAULT/DECLARED GLOBAL VARIABLES
######################

my $deldup_window = 1000000;  # default 1m
my $sv_close="100000";

my $pid_list;
my @pids;
my $file_list;
my $output_prefix;
my $help;

######################
## PARSE INPUTS
######################

GetOptions("results_dir=s"         => \$results_dir,
           "sv_close=i"            => \$sv_close,
           "deldup_window=i"       => \$deldup_window,
           "pids=s"                => \$pid_list,
           "file_list=s"           => \$file_list,
	   "output_prefix=s"       => \$output_prefix,
           "help"                  => \$help
           ) or die ("Error in command line arguments\n");

if ($help){
  die ($usage);
}

warn "\n";
warn "WARNING: output_prefix not defined\n"                       unless defined $output_prefix;
warn "WARNING: could not find results directory: $results_dir\n"  unless -e $results_dir;
warn "WARNING: could not find genome fai file: $genome_fai\n"     unless -e $genome_fai;
warn "WARNING: could not find genes bed file: $genes_bed\n"       unless -e $genes_bed;

die ($usage) unless defined $output_prefix;
die ($usage) unless -e $results_dir;
die ($usage) unless -e $genome_fai;
die ($usage) unless -e $chrom_arms;
die ($usage) unless -e $genes_bed;

if (length $pid_list > 1){
  warn "WARNING: variable \$pid_list was defined as \"$pid_list\"...\n";
  if (-e $pid_list){
       warn "WARNING: variable \$pid_list was defined as \"$pid_list\"... which is a file... using PID definition in file\n";
       open(my $pid_fh, "<$pid_list") or die "ERROR: cannot open file $pid_list";	
       while(my $pid_line = <$pid_fh>){
         chomp $pid_line;
         my @pids_line_arr = split (" ",$pid_line);
         push(@pids, @pids_line_arr);
       }
       close($pid_fh);
       warn "\#\# PIDs found in $pid_list: ".join(" ",@pids)."\n";
  }
  else {
       warn "WARNING: variable \$pid_list was defined as \"$pid_list\"... which is not a file... parsing space separated PIDs\n";
       @pids = split (" ",$pid_list);
  }
}
else{
  if (length $file_list > 1){
    warn "WARNING: file list defined... using \"$file_list\"\n";
     open(my $pid_fh, "<$file_list") or die "ERROR: cannot open file $file_list";
     while(my $pid_line = <$pid_fh>){
       chomp $pid_line;
       push(@pids, $pid_line);
     }
     close($pid_fh);
     warn "\#\# PIDs found in $pid_list: ".join(" ",@pids)."\n";

  }
  else{
    warn "WARNING: no pid list supplied... using all pid is results folder \"$results_dir\"\n";
    @pids = `ls $results_dir`
  }
}

chomp @pids;
my $num_pids = scalar @pids;

my $outfile             = "$output_prefix.$version.svNear$sv_close.deldupwin$deldup_window.SV.for_circlize.tsv";
my $files_used          = "$output_prefix.$version.svNear$sv_close.deldupwin$deldup_window.SV.files_used.tsv";
my $log_file            = "$output_prefix.$version.svNear$sv_close.deldupwin$deldup_window.SV.log.txt";

##########################
## print the parameters ##
##########################

warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";
warn "\# PROJECT DIR AND PIDS \#\n";
warn "GENOME_FAI:     $genome_fai\n";
warn "GENES_BED:      $genes_bed\n";
warn "RESULTS_DIR:    $results_dir\n";
warn "PIDS:           @pids\n";
warn "NUM_PIDS:       $num_pids\n";
warn "\n";
warn "\# FOLDER/FILE WILDCARDS \#\n";
warn "SV FILE:        PID/$folder_sv/$pre_sv"."[PID]"."$ext_sv\n";
warn "\n";
warn "\# VARIABLES \#\n";
warn "SV_CLOSE:       $sv_close\n";
warn "DELDUP_WINDOW:  $deldup_window\n";
warn "\n";
warn "\# OUTFILES \#\n";
warn "GENERECURRENCE: $outfile\n";
warn "FILES USED:     $files_used\n";
warn "LOG FILE:       $log_file\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

######################
## Check for SV files
######################

my %mutation_files;
my @sv_file_list;

warn "\# CHECKING: PIDS, SAMPLES and SV files ...\n";

my @pid_samples_temp;

open (my $files_used_fh, ">$files_used") or die "\n\nERROR: cannot open file: $files_used\n\n";
open (my $log_file_fh, ">$log_file")     or die "\n\nERROR: cannot open file: $log_file\n\n";

### Looping over the pids

if (length $file_list > 1){
  foreach my $file (@pids){
    chomp $file;
    $file =~ m/^.*?svs_(.+?)_filtered_somatic_minEventScore3.tsv.*?/;
    my $pid = $1;
    my $pid_sample = $pid;
    my $sample_combination = $pid;

    print $log_file_fh "\nProcessing: $pid\t$sample_combination\n";

    if (-e $file){
      warn "Processing: $pid\t[FILE: \"$file\"]\n";
      $mutation_files{$pid_sample}{"pid"}                    = "$pid";
      $mutation_files{$pid_sample}{"sample"}                 = "$pid";
      $mutation_files{$pid_sample}{"sv_file"}                = "$file";
      print $files_used_fh "$pid\t$sample_combination\t$pid_sample\tSV\t".$mutation_files{$pid_sample}{"sv_file"}."\n";
      push (@pid_samples_temp, $pid_sample);
    }
    else{
      print $log_file_fh "-- WARNING: could not find SV file: ".$file."\n";
      warn "Processing: $pid\t[FILE: \"$file\"] [no SV]\n";
    }
  }
}
else{
foreach my $pid (@pids){
  chomp $pid;

  my @SV_dirs = `ls $results_dir/$pid/$folder_sv* 2>/dev/null`;
  if (length($SV_dirs[0]) < 3 ){
    warn "WARNING: cannot find SV dir for $pid ($results_dir/$pid/$folder_sv*\*) .... skipping\n";
    next;
  }
  
  chomp(my @sample_combinations = `ls -d $results_dir/$pid/$folder_sv* | xargs -n1 basename`);

  my @sample_combinations_new = map {s/$folder_sv.*?//r; } @sample_combinations;

  foreach my $sample_combination(@sample_combinations_new) {

    print $log_file_fh "\nProcessing: $pid\t$sample_combination\n";
    my $pid_sample = $pid."_".$sample_combination;
    push (@pid_samples_temp, $pid_sample);

    my %expected_files;
    $expected_files{"sv_file"}="$results_dir/$pid/$folder_sv$sample_combination/$pre_sv$pid$ext_sv";

    if (-e $expected_files{"sv_file"}){
      warn "Processing: $pid\t$sample_combination\n";
      $mutation_files{$pid_sample}{"pid"}                    = "$pid";
      $mutation_files{$pid_sample}{"sample"}                 = "$sample_combination";
      $mutation_files{$pid_sample}{"sv_file"}                = (`ls $expected_files{"sv_file"} 2>/dev/null`);
      chomp $mutation_files{$pid_sample}{"sv_file"};
      print $files_used_fh "$pid\t$sample_combination\t$pid_sample\tSV\t".$mutation_files{$pid_sample}{"sv_file"}."\n";
    }
    else{
      warn "Processing: $pid\t$sample_combination [no SV]\n";
      print $log_file_fh "-- WARNING: could not find SV file: ".$expected_files{"sv_file"}."\n";
    }
  }
}
}

close ($files_used_fh);

warn "\n\#\# Found ".scalar(@pid_samples_temp)." samples from ".scalar(@pids)." pids... switching pid-sample identifiers...\n\n";
print $log_file_fh "\n\n-- Found ".scalar(@pid_samples_temp)." samples from ".scalar(@pids)." pids... switching pid-sample identifiers...\n";

@pids = @pid_samples_temp;

warn "\# DONE! Found all files\n";
print $log_file_fh  "-- DONE! Found all files\n\n";

close($log_file_fh);

######################
## OUTPUT DATASTRUCTURE DEFINITIONS
######################

my $sample_number = 1;
my %window_mapping=();
my %window_genes=();
my %gene_mapping=();
my %gene_mapping_exact=();
my %gene_mapping_close=();
my %del_mapping=();
my %dup_mapping=();
# $sv_close - varaible containing the threshold of a "close" SV

######################
## RUN DETAILS
######################

warn "\n";
warn "# RUN :: $0 @arguments\n";
warn "# RUN :: num_files = " . scalar(@pids) . "\n";

######################
## LOOP OVER FILES
######################

foreach my $pid_sample (@pids){

  my $sv_file =  $mutation_files{$pid_sample}{"sv_file"};

  #########################
  # Basic sample statistics
  ##########################

  warn "\n# SAMPLE $sample_number :: file = $sv_file\n";
  my $sv_file_lines = ` grep -v ^# $sv_file | wc -l | cut -f 1 -d " "`;
  chomp $sv_file_lines;
  warn "# SAMPLE $sample_number :: num_SVs = $sv_file_lines\n";
  my $sv_file_TRA_lines = ` grep -v ^# $sv_file | awk '\$9=="TRA"' | wc -l | cut -f 1 -d " "`;
  my $sv_file_INV_lines = ` grep -v ^# $sv_file | awk '\$9=="INV"' | wc -l | cut -f 1 -d " "`;
  my $sv_file_DUP_lines = ` grep -v ^# $sv_file | awk '\$9=="DUP"' | wc -l | cut -f 1 -d " "`;
  my $sv_file_DEL_lines = ` grep -v ^# $sv_file | awk '\$9=="DEL"' | wc -l | cut -f 1 -d " "`;
  chomp $sv_file_TRA_lines;
  chomp $sv_file_INV_lines;
  chomp $sv_file_DEL_lines;
  chomp $sv_file_DUP_lines;
  my $sv_file_TRAINVDUPDEL_lines = $sv_file_TRA_lines + $sv_file_INV_lines + $sv_file_DUP_lines + $sv_file_DEL_lines;
  warn "# SAMPLE $sample_number :: num(DEL,DUP,INV,TRA,total) = ($sv_file_DEL_lines,$sv_file_DUP_lines,$sv_file_INV_lines,$sv_file_TRA_lines,$sv_file_TRAINVDUPDEL_lines)\n";

  if ($sv_file_TRAINVDUPDEL_lines ne $sv_file_lines){ 
    die "# ERROR :: SAMPLE $sample_number :: number of SVs not consistent ($sv_file_TRAINVDUPDEL_lines ne $sv_file_lines)";
  }

  #########################
  # STORE SVs
  #########################

  my @sv_part1 = `grep -v ^# $sv_file | cut -f 1-3`;
  my @sv_part2 = `grep -v ^# $sv_file | cut -f 4-6`;

  foreach my $sv1 (@sv_part1){
    my $sv2 = shift(@sv_part2);

    chomp $sv1;
    chomp $sv2;

    $window_mapping{$sv1}{$sv2}{$sample_number} = 1;
  }  

  #########################
  # Assign genes to SVs per sample
  #########################

  # create a bedtools genome file... an ugly solution but it works

  my $temp_genome = "___genome.txt";
  `cut -f 1,2 $genome_fai >  $temp_genome`;

  # identify closest gene (currently protein coding genes) to each SV using closestbed

  my @close_gene_list_1 = `grep -v ^# $sv_file | cut -f 1-3 | bedtools sort -faidx $genome_fai -i - | closestBed -g $temp_genome -d -a - -b $genes_bed | cut -f 1-7,13`;
  my @close_gene_list_2 = `grep -v ^# $sv_file | cut -f 4-6 | bedtools sort -faidx $genome_fai -i - | closestBed -g $temp_genome -d -a - -b $genes_bed | cut -f 1-7,13`;
  `rm  $temp_genome`;

  foreach my $gene1 (@close_gene_list_1){
    my $gene2 = shift(@close_gene_list_2);

    chomp $gene1;
    chomp $gene2;

    my @rray_gene1 = split("\t", $gene1);
    my @rray_gene2 = split("\t", $gene2);

    my $pos1 = join ("\t", $rray_gene1[0], $rray_gene1[1], $rray_gene1[2]);
    my $pos2 = join ("\t", $rray_gene2[0], $rray_gene2[1], $rray_gene2[2]);

    my $gene_pos1 = join ("\t", $rray_gene1[3], $rray_gene1[4], $rray_gene1[5], $rray_gene1[6]);
    my $gene_pos2 = join ("\t", $rray_gene2[3], $rray_gene2[4], $rray_gene2[5], $rray_gene2[6]);

    $gene_mapping{$gene_pos1}{$sample_number}=1;
    $gene_mapping{$gene_pos2}{$sample_number}=1;
    
    $gene_mapping_exact{$gene_pos1}{$sample_number}=1 if $rray_gene1[7] < 1;
    $gene_mapping_exact{$gene_pos2}{$sample_number}=1 if $rray_gene2[7] < 1;

    $gene_mapping_close{$gene_pos1}{$sample_number}=1 if $rray_gene1[7] <= $sv_close;
    $gene_mapping_close{$gene_pos2}{$sample_number}=1 if $rray_gene2[7] <= $sv_close;

    $window_genes{$pos1}{"gene1"} = $gene_pos1;
    $window_genes{$pos2}{"gene2"} = $gene_pos2;
  }

  #########################
  # identify deletion regions
  #########################

  my @deletion_list = `grep -v ^# $sv_file | awk -F'\t' '\$9 == "DEL"' | cut -f 1,2,6`;
  foreach my $deletion (@deletion_list){
    chomp $deletion;
    my @rray_deletion = split ("\t", $deletion);
    if ($rray_deletion[1] > $rray_deletion[2]){
      ($rray_deletion[1],$rray_deletion[2]) = ($rray_deletion[2],$rray_deletion[1]);
    }
    $rray_deletion[1] = (int($rray_deletion[1]/$deldup_window))*$deldup_window + 1;
    $rray_deletion[2] = (int($rray_deletion[2]/$deldup_window))*$deldup_window + $deldup_window;

    while ($rray_deletion[1] < $rray_deletion[2]){  
      my $end = ($rray_deletion[1] + $deldup_window - 1);
      $del_mapping{"$rray_deletion[0]\t$rray_deletion[1]\t$end"}{$sample_number} = $del_mapping{"$rray_deletion[0]\t$rray_deletion[1]\t$end"}{$sample_number} + 1;
      $rray_deletion[1] = $rray_deletion[1] + $deldup_window;
    }
  }

  #########################
  # identify duplication regions
  #########################

  my @duplication_list = `grep -v ^# $sv_file | awk -F'\t' '\$9 == "DUP"' | cut -f 1,2,6`;
  foreach my $duplication (@duplication_list){
    chomp $duplication;
    my @rray_duplication = split ("\t", $duplication);
    if ($rray_duplication[1] > $rray_duplication[2]){
      ($rray_duplication[1],$rray_duplication[2]) = ($rray_duplication[2],$rray_duplication[1]);
    }
    $rray_duplication[1] = (int($rray_duplication[1]/$deldup_window))*$deldup_window + 1;
    $rray_duplication[2] = (int($rray_duplication[2]/$deldup_window))*$deldup_window + $deldup_window;

    $dup_mapping{"$rray_duplication[0]\t$rray_duplication[1]\t$rray_duplication[2]"}{$sample_number} = $dup_mapping{"$rray_duplication[0]\t$rray_duplication[1]\t$rray_duplication[2]"}{$sample_number} + 1;
  }

  #########################
  # end loop, incrament counter
  #########################

  $sample_number = $sample_number + 1;
}

#########################
# PRINT HEADER
#########################

open (my $outfile_fh, ">$outfile") or die;
print $outfile_fh "TYPE\tCHR_1\tSTART_1\tEND_1\tCHR_2\tSTART_2\tEND_2\tNAME\tSCORE\tSTRAND\tCOMMENT\n";

#########################
# PRINT INV and TRA SVs in windows
#########################

## print: position 1, position 2, uniqe samples with SVs, total SV count

foreach my $pos1 (sort keys %window_mapping){
  chomp $pos1;
  foreach my $pos2 (sort keys %{ $window_mapping{$pos1} }){
    chomp $pos2;

    my $gene1_pos =  $window_genes{$pos1}{"gene1"};
    my $gene2_pos =  $window_genes{$pos2}{"gene2"};

    chomp $gene1_pos;
    chomp $gene2_pos;
   
    my $gene1_score = scalar(keys %{ $gene_mapping{$gene1_pos} });
    my $gene2_score = scalar(keys %{ $gene_mapping{$gene1_pos} });
    my $gene1_score_exact = scalar(keys %{ $gene_mapping_exact{$gene1_pos} });
    my $gene2_score_exact = scalar(keys %{ $gene_mapping_exact{$gene1_pos} });
    my $gene1_score_close = scalar(keys %{ $gene_mapping_close{$gene1_pos} });
    my $gene2_score_close = scalar(keys %{ $gene_mapping_close{$gene1_pos} });    

    my $max_score       = $gene1_score;
    my $max_score_exact = $gene1_score_exact;
    my $max_score_close = $gene1_score_close;

    $max_score = $gene2_score if ($gene2_score > $gene1_score);
    $max_score_exact = $gene2_score_exact if ($gene2_score_exact > $gene1_score_exact);
    $max_score_close = $gene2_score_close if ($gene2_score_close > $gene1_score_close);
   
    print $outfile_fh "windows_SVs_direct\tchr$pos1\tchr$pos2\t.\t$max_score_exact\t.\t.\n" unless (length ($pos2) < 3);
    print $outfile_fh "windows_SVs_near\tchr$pos1\tchr$pos2\t.\t$max_score_close\t.\t.\n" unless (length ($pos2) < 3);
    print $outfile_fh "windows_SVs_closest\tchr$pos1\tchr$pos2\t.\t$max_score\t.\t.\n" unless (length ($pos2) < 3);
  }
}

#########################
# PRINT SVs in genes
#########################

## print: position 1, gene name, unniqe samples with SVs, total SV count for gene

foreach my $gene (sort keys %gene_mapping){
  chomp $gene;
  my $gene_sample_count = scalar(keys %{ $gene_mapping{$gene}}) ;
  my $gene_sv_count = 0;
  foreach my $sample_sv (keys %{ $gene_mapping{$gene}}){
    $gene_sv_count = $gene_sv_count + $gene_mapping{$gene}{$sample_sv};
  }
  my @gene_split = split("\t", $gene);
  print $outfile_fh "gene_SVs_closest\tchr$gene_split[0]\t$gene_split[1]\t$gene_split[2]\t.\t0\t0\t$gene_split[3]\t$gene_sample_count\t.\t$gene_split[3] ($gene_sample_count)\n" unless (length ($gene_split[0]) < 1) ;
}

foreach my $gene (sort keys %gene_mapping_exact){
  chomp $gene;
  my $gene_sample_count = scalar(keys %{ $gene_mapping_exact{$gene}}) ;
  my $gene_sv_count = 0;
  foreach my $sample_sv (keys %{ $gene_mapping_exact{$gene}}){
    $gene_sv_count = $gene_sv_count + $gene_mapping_exact{$gene}{$sample_sv};
  }
  my @gene_split = split("\t", $gene);
  print $outfile_fh "gene_SVs_direct\tchr$gene_split[0]\t$gene_split[1]\t$gene_split[2]\t.\t0\t0\t$gene_split[3]\t$gene_sample_count\t.\t$gene_split[3] ($gene_sample_count)\n" unless (length ($gene_split[0]) < 1) ;
}

foreach my $gene (sort keys %gene_mapping_close){
  chomp $gene;
  my $gene_sample_count = scalar(keys %{ $gene_mapping_close{$gene}}) ;
  my $gene_sv_count = 0;
  foreach my $sample_sv (keys %{ $gene_mapping_close{$gene}}){
    $gene_sv_count = $gene_sv_count + $gene_mapping_close{$gene}{$sample_sv};
  }
  my @gene_split = split("\t", $gene);
  print $outfile_fh "gene_SVs_near\tchr$gene_split[0]\t$gene_split[1]\t$gene_split[2]\t.\t0\t0\t$gene_split[3]\t$gene_sample_count\t.\t$gene_split[3] ($gene_sample_count)\n" unless (length ($gene_split[0]) < 1) ;
}


#########################
# PRINT DELs in dels
#########################

## print: position 1, del name, unniqe samples with SVs, total SV count for del

foreach my $del (sort keys %del_mapping){
  chomp $del;
  my $del_sample_count = scalar(keys %{ $del_mapping{$del}});
  my $del_sv_count = 0;
  foreach my $sample_sv (keys %{ $del_mapping{$del}}){
    $del_sv_count = $del_sv_count + $del_mapping{$del}{$sample_sv};
  }
  print $outfile_fh "del_SVs\tchr$del\t.\t0\t0\t.\t$del_sample_count\t.\t.\n" ;
}

#########################
# PRINT DUPs in dups
#########################

## print: position 1, dup name, unniqe samples with SVs, total SV count for dup

foreach my $dup (sort keys %dup_mapping){
  chomp $dup;
  my $dup_sample_count = scalar(keys %{ $dup_mapping{$dup}});
  my $dup_sv_count = 0;
  foreach my $sample_sv (keys %{ $dup_mapping{$dup}}){
    $dup_sv_count = $dup_sv_count + $dup_mapping{$dup}{$sample_sv};
  }
  print $outfile_fh "dup_SVs\tchr$dup\t.\t0\t0\t.\t$dup_sample_count\t.\t.\n" ;
}

close($outfile_fh);

warn "\nCompleted successfully!\n";
