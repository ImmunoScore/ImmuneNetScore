#!/usr/bin/perl

use warnings;
use IO::File;
use Cwd;
use lib "lib" ;
use ImmuneScore;
use ImmuneGEO;
use ImmuneCluster;

# quit unless correct number of command-line args 
my $num_args = $#ARGV + 1;
if ($num_args != 3) {
  print "\nUsage: ImmuneProfiler.pl, lacking commmand line arguments\n";
  print "arg1: Case name\n";
  print "arg2: Quantile normalized GEO matrix file name\n";
  print "arg3: Configuration file for the GEO matrix case to run ImmuneScore\n";
  exit;
}
my ($case, $file_config, $file_GEO) =@ARGV;

# opend file handles
my $fh_config = new IO::File "$file_config" or die "$file_config $!";
my $fh = new IO::File "$file_GEO" or die "$file_GEO  $!";

# configuration file and properties
my $cwd = getcwd();
my $DATADIR = "$cwd/data";
print "Processing the config file for $case case. Data directory: $DATADIR \n";

# fetch config file data
my ($rh_config,$rh_properties) = ImmuneGEO::config_immune($fh_config);
my $GEO_source = $rh_config->{accession}; 
my $annotation_name = $rh_config->{platform};
my $gene_col_platform = $rh_config->{gene_platform};
my $begin_platforn = $rh_config->{begin_data_platform};
my $patient_row = $rh_config->{patient_id_row};
my $sample_row = $rh_config->{sample_id_row};
my $exclude = $rh_config->{exclude}; 
my $saturation_limit = $rh_config->{saturation} || 0.9;
my $logratio_limit = $rh_config->{logratio} || 0.5; 

#configure files and folders
my $saturation_DIR = "$DATADIR\/saturation";
unless(-e $saturation_DIR ) {
   die "Unable to open $saturation_DIR\n";
}

my $apcuster_DIR = "$DATADIR\/apclusters";
unless(-e $apcuster_DIR ) {
   die "Unable to open $apcuster_DIR\n";
}

my $platform_DIR = "$DATADIR\/platform";
unless(-e $platform_DIR ) {
   die "Unable to open $platform_DIR\n";
}

my $case_DIR = "$DATADIR\/$case";
unless(-e $case_DIR or mkdir $case_DIR) {
   die "Unable to create $case_DIR\n";
}

my $result_dir = "$case_DIR\/results";
unless(-e $result_dir or mkdir $result_dir) {
   die "Unable to create $result_dir\n";
}

#file declarations and file handles
my $fh_names = new IO::File "$DATADIR\/Immune_cell_names.txt" or die "cannot open $DATADIR\/Immune_cell_names.txt $!";
my $fh_names_res = new IO::File "> $DATADIR\/attribute_names.txt" or die "cannot open $DATADIR\/attribute_names.txt $!";
my $fh_net = new IO::File "$DATADIR\/immune_cell_network.txt" or die "cannot open $DATADIR\/immune_cell_network.txt $!";
my $fh_net_res = new IO::File "> $DATADIR\/immune_network.txt" or die "cannot open $DATADIR\/immune_network.txt $!";
my $file_marker_genes = "$case_DIR\/DIST_Marker_genes_$case\_subtype_scores\.txt";
my $fh_marker_genes = new IO::File "> $file_marker_genes" or die "$file_marker_genes does not exist $!";
my $file_marker_genes_GIT = "$case_DIR\/GIT_Marker_genes_$case\_subtype_scores.txt";
my $fh_marker_genes_GIT = new IO::File "> $file_marker_genes_GIT" or die "$file_marker_genes_GIT does not exist\n";

# declarations
my @subtypes = qw(Macrophage Dendritic cytotoxicT Th1 Th2 treg natural_killer B_cells Granulocytes MDSC Mast_cells); # GIT 
my %patient2subtypescore = () ;
my %subtype2patient_data=();
my @row_ids = ();
my @headers = ();
my %matrix = ();
my %columns_data = ();
my $row_index = 0;
my %header_column_index =();
my $data_col_end;
my %gene_rank_count=();
my %gene_expression_sum=();
my %code2name = ();
my %immune_cell_net = ();
my $ra_expanded_subtypes;
my $rh_exclude;
my @expanded_subtypes = ();
my @network_informed_subtypes = ();
my @general_subtypes = ();
my %total_subtype_expression = (); 

# immune cell cytokine network code name mapping
my %immune_name2code = ('cytotoxicT'=>85,
                        'natural_killer'=>8,
                        'Th1'=>83,
                        'Th2'=>84,
                        'B_cells'=>17,
                        'Dendritic'=>146,
                        'treg'=>191,
                        'Macrophage'=>2,
                        'Mast_cells'=>88,
                        );

print "Calculating the General (GIT) and distinct (DIST) immune-cell subset scores, and immune-cell networks...\n";
# fetch the immunce cell phenotype speciifc ap clusters calculated from the PPI similarity network of the cells
my ($rh_specific_type2gene_apclusters,$rh_specific_type2gene_cluster_IDs) = ImmuneCluster::specific_type_apclusters($apcuster_DIR);
 

# fetch the semantic literature saturation point for the genes
my ($rh_saturation, $rh_subtype2genes) = ImmuneScore::fetch_saturation(\@subtypes,
                                                                        $saturation_DIR,
                                                                        $fh_marker_genes_GIT,
                                                                        $saturation_limit
                                                                        );

# network informed saturation scores and expanded subtypes
my $rh_cell_specific_networks = ImmuneScore::cell_specific_networks();
my %cell_specific_networks = %$rh_cell_specific_networks;

# immune networks node attributes
print $fh_names_res "code,name\n"; 
while (! $fh_names->eof()){
  my $line;
  while (defined ($line = <$fh_names>)){
    $line  =~s/\s+$//; $line  =~s/"//g;
    my @line = split("\,",$line);
    my $code = $line[0];  $code  =~s/\s+$//; $code  =~s/"//g;
    my $name = $line[1];
    $code2name{$code} = $name; 
    print $fh_names_res "$code,$name\n"; 
  }
}

# build GIT template network
while (! $fh_net->eof()){
  my $line;
  while (defined ($line = <$fh_net>)){
    $line  =~s/\s+$//; $line  =~s/"//g;
    $line  =~s/\s+/,/g;
    my @line = split("\,",$line);
    my $n1 = $line[0];$n1  =~s/\s+$//; $n1  =~s/"//g;
    my $n2 = $line[1];$n2  =~s/\s+$//; $n2  =~s/"//g;
    next unless (exists $code2name{$n1} && exists $code2name{$n2});
    my $name1= $code2name{$n1};
    my $name2= $code2name{$n2};
    $immune_cell_net{$n1}{$n2}++;
    $immune_cell_net{$n2}{$n1}++; 
    print $fh_net_res "$n1,$n2\n"; 
  }
}

# Expand the immune cell marler list with Protein-protein interaction informed Affinity Purification clusters
($rh_saturation, $ra_expanded_subtypes,$rh_exclude,$rh_DIST_to_GIT) = ImmuneScore::network_informed_saturation($rh_specific_type2gene_apclusters,
                                                                                                               $rh_specific_type2gene_cluster_IDs,
                                                                                                               $rh_saturation,
                                                                                                               \@subtypes,
                                                                                                               \%cell_specific_networks,
                                                                                                               $rh_subtype2genes,
                                                                                                               $fh_marker_genes
                                                                                                              );

foreach my$type(@$ra_expanded_subtypes){
  next if exists $rh_exclude->{$type};
  push (@expanded_subtypes,$type);
  if ($type  =~ /^EXT/){
       push (@network_informed_subtypes,$type)  
  }else{
      push (@general_subtypes,$type)  
  } 
}
@subtypes = @expanded_subtypes;
my %stratified_subtypes = ('general'=>\@general_subtypes,
                           'specific'=> \@network_informed_subtypes);

# process gene annotations
my $rh_illumina2gene = ImmuneGEO::fetch_platform_annotations($platform_DIR,
                                                             $annotation_name,
                                                             $begin_platforn, 
                                                             $gene_col_platform 
                                                             )  ; 


# Building data matrix from the quintile normalized GEO file...\n";
while (! $fh->eof()){
  my $line;
  while (defined ($line = <$fh>)){
    $line  =~s/\s+$//; $line  =~s/"//g;
    $i++;
    1; 
    if ($i==1){
        @header_row = split("\t", $line);
        my $h = 0;
        my $k=0; 
        foreach my$header(@header_row){
          $k++;
          if ($header =~/GSM\d+/) {
             $data_col_end = $k ;
             $h++;
             push (@headers,$header) ;
             $header_column_index{$h} = $header;
          } 
        } 
        next; 
    } 
  
    my @line= split("\t", $line); 
    my $id = $line[0];  $id  =~s/\s+$//;
    push (@row_ids, $id);
    my $j=0;
    my $column_index = 0 ; 
    $row_index++; 
    foreach my$feature(@line){
       $j++;
       next if $j ==1;
       last if $line =~ /^(\!)/;
       $column_index++ ;
       ($feature eq 'NaN' || $feature eq 'NA'  ) ? $feature = 0 : $feature = $feature;
       $matrix{$row_index}{$column_index} = $feature;
       push(@{$columns_data{$column_index}},[$feature,$row_index]);
    } 
  } 
}

#gene expression profiles for each patient/sample 
my %row_data=();
foreach my$header(@headers){ 
  $n++;
  my @columndata = @{$columns_data{$n}};
  my $s = scalar (@columndata);
  my @row_index = map {$_->[1]} @columndata;
  my @features = map {$_->[0]} @columndata;

  my $r = 0;
  foreach my$row_index(@row_index){
    $r++;
    my $feature = $matrix{$row_index}{$n} ;
    my $id = $row_ids[$row_index-1];
    my $sym = $rh_illumina2gene->{$id} || '--'; $sym  =~s/\s+$//;
    push(@{$row_data{$sym}},$feature); 
  } 

  my $rank=0;
  my %genes = ();
  foreach my$row_index(@row_index){
    $rank++;
    my $feature = $matrix{$row_index}{$n} || 0;
    my $id = $row_ids[$row_index-1];
    my $sym = $rh_illumina2gene->{$id} || '--';
    
    # retreive the immune subtypes for the gene
    foreach my$subtype(@subtypes){
      my $ra_score =  $rh_saturation->{$sym}{$subtype} ;
      my $info = $ra_score->[0];
      my $saturation = $ra_score->[1] || 0;
    
      # select gene expression rank which approaches a saturation 
      if ((defined($saturation)) && ($saturation > $saturation_limit)){
        # count the number of genes contributing to the sbutype ranksum for each patient/sample
        $gene_rank_count{$subtype}{$n}++ ; 
        
        # sum the total  number of gene expression contributing to the sbutype for each paitent/column n
        my $previous_exp =  $gene_expression_sum{$subtype}{$n} || 0; 
        $gene_expression_sum{$subtype}{$n} = $previous_exp + $feature;
        my $new_exp =  $gene_expression_sum{$subtype}{$n} ;
        
        # populate a data structure which stores the subtype marker genes for each subtype across patiens
        $feature eq '-inf' ? $feature = 0 : $feature = $feature;
        
      }  
    }# end subtype 
  } # end gene 
} # end patient column / header

# the average gene expression value for each subtype,over all patients/columns...\n";
my $p=0;
foreach my$header(@headers){
  $p++;
  foreach my $subtype(@subtypes){
    # look up the ranksum for the subtype for patient "p"
    my $expression_sum = $gene_expression_sum{$subtype}{$p}; 

    # calculate the average expression over all genes contributing to the subtype total expression score
    my $gene_rank_count = $gene_rank_count{$subtype}{$p} || 0;
    my $avg_expression ;
    $gene_rank_count == 0 ?   $avg_expression = 0  : $avg_expression = ($expression_sum/$gene_rank_count);  

    #calculate the total avg for the subtype across all patients
    my $previous_expr = $total_subtype_expression{$subtype} ||0 ;
    $total_subtype_expression{$subtype} = $previous_expr + $avg_expression;
  }#end foreach subype
   
} # end foreach header columns/patient


# stratify for general and specific subtypes (GIT or DIST )s
foreach my$type(keys %stratified_subtypes){
   my @subtypes = @{$stratified_subtypes{$type}};
   my $res = "$type\_res\_$case\_subtype_scores\.txt";
   my $res_FILE = "$result_dir\/$res";
   my $fh_res = new IO::File "> $res_FILE" or die "cannot open $res_FILE $!";
  
   my $subytpes = join("\t",@subtypes);$subytpes =~s/EXT_//g; 
   print $fh_res "GEO accession \t$subytpes\n";
   
   my $q = 0;
   foreach my$header(@headers){
    $q++;
    my @patient_norm_subtypes = ();
    foreach my $subtype(@subtypes){
     
      my $expression_sum = $gene_expression_sum{$subtype}{$q}; 
      # total exoression signature across the subtype
      my $total_expr =  $total_subtype_expression{$subtype} || 0; 
      my $n_patients = scalar(@headers);
     
      # calculate the average expression over all genes contributing to the subtype total expression score
      my $gene_rank_count = $gene_rank_count{$subtype}{$p} || 0;
      my $avg_expression ;
      $gene_rank_count == 0 ?   $avg_expression = 0  : $avg_expression = ($expression_sum/$gene_rank_count);  
  
      #avg total expression for the subtype  (over all patients)
      my $avg_patient_expression = $total_expr / $n_patients;
      # normalized expression (normalized expression for the subtype signature genes, normalized for the average per patient)
      my $norm_expression ;
      $avg_patient_expression == 0 ?   $norm_expression = 0  :  $norm_expression = ($avg_expression/$avg_patient_expression);
      ($norm_expression == -0.0 || $norm_expression eq '-nan' ) ?     $norm_expression = 0 : $norm_expression = $norm_expression;
      
      push (@patient_norm_subtypes,$norm_expression);
      $patient2subtypescore{$header}{$subtype}= [$norm_expression,$q];
      push(@{$subtype2patient_data{$subtype}}, [$norm_expression,$header,$q]);
    }#end foreach subype
    
    my $patient_subytpes = join("\t",@patient_norm_subtypes);
    print $fh_res "$header\t$patient_subytpes\n" ;
  } # end foreach header columns/patient
   
  $fh_res -> close;
  $fh -> close;
  
  # Ratio scores for the immune subtypes 
  print  "compute the subtype ratios(subtype: $type)\n" ;
  my ($rh_patient2ratio_score, $ra_subtype_ratios,$rh_subtype_ratio2scores) = ImmuneScore::subset_ratio_network(\%patient2subtypescore,
                                                                                                                \@headers,
                                                                                                                \@subtypes,
                                                                                                                $type,
                                                                                                                $result_dir,
                                                                                                                $rh_DIST_to_GIT,
                                                                                                                \%code2name,
                                                                                                                \%immune_name2code,
                                                                                                                \%immune_cell_net,
                                                                                                                $case,
                                                                                                                $logratio_limit);
   
} # end fore each stratified immne subtype groups


print "DONE: Immune-cell subset scores and networks generated\n";
exit ;



