##############################################################################
# file: ImmuneScore.pm
# desc: ImmuneScore package
##############################################################################

package ImmuneScore;

use strict;
use IO::File;
use List::Util;
use POSIX qw(ceil floor);
 
sub network_informed_saturation {

 my $rh_subtype2gene_apclusters = shift;
 my $rh_subtype2gene_cluster_IDs = shift;
 my $rh_saturation = shift;
 my $ra_subtype = shift;
 my $rh_cell_specific_networks = shift;
 my $rh_subtype2genes = shift;
 my $fh_marker_genes = shift;
 
 my %subtype2gene_apclusters = %$rh_subtype2gene_apclusters;
 # map of the DIST to the GIT
 my %specific_type2subtype = ();
 my %DIST_to_GIT = ();
 while(my($subtype, $ra_specific_type)=each(%$rh_cell_specific_networks)){
    foreach(@$ra_specific_type){ 
       $specific_type2subtype{$_} =  $subtype;
    }
 }
 
 my %saturation = %$rh_saturation;
 my @subtypes = @$ra_subtype;
 my %exclude = ();
 foreach my $specific_type(keys %subtype2gene_apclusters){
    next unless exists $specific_type2subtype{$specific_type};
    my $subtype = $specific_type2subtype{$specific_type};
    
    # assign a new subtype name which maps the subtype to the specfifc type
    my $expanded_subtype = "EXT_$subtype\_$specific_type";
    $expanded_subtype  =~s/\s+/_/g;
    $expanded_subtype  =~s/\+/Pos/g;
    $expanded_subtype  =~s/\-/Neg/g;
    
    $DIST_to_GIT{$expanded_subtype} = $subtype; 
    
    my @subtype_genes = @{$rh_subtype2genes->{$subtype}};
    foreach my$subt_gene(@subtype_genes){
       my $saturation = $rh_saturation->{$subt_gene}{$subtype}->[1];
       my $info =  $rh_saturation->{$subt_gene}{$subtype}->[0] || 0;
       my $cocites =  $rh_saturation->{$subt_gene}{$subtype}->[2] || 0;
       $saturation{$subt_gene}{$expanded_subtype} = [$info,$saturation,$cocites] ;
   }

    push(@subtypes,$expanded_subtype) ;
    
    #flag to exclude specific cell type if they do not have an ap cluster
    my $exclude = 1;
    
    my $ra_clusters = $rh_subtype2gene_cluster_IDs->{$specific_type};
    my %saturated_in_cluster = ();
    foreach my$cluster_ID(@$ra_clusters){ 
       my $rh_gene2clustergenes = $subtype2gene_apclusters{$specific_type}{$cluster_ID};
    
       foreach my$gene(keys %$rh_gene2clustergenes){
          my $saturation = $rh_saturation->{$gene}{$subtype}->[1];
          my $info =  $rh_saturation->{$gene}{$subtype}->[0] || 0;
          my $cocites =  $rh_saturation->{$gene}{$subtype}->[2] || 0;
          
          my $ra_genes = $rh_gene2clustergenes->{$gene};
          if ($saturation >= 0.9){
             $exclude = 0; 
             $saturated_in_cluster{$gene}++; 
             $saturation{$gene}{$expanded_subtype} = [$info,$saturation,$cocites] ;
             print $fh_marker_genes "$gene\t$gene\t$specific_type\t$subtype\tGIT\n" ;
             foreach my$gene_ap(@$ra_genes){
               my $type = '';
               if (exists($saturated_in_cluster{$gene_ap})) {
                  $type = 'GIT';
               }else{
                  $type = 'DIST';
               }

               my $sat =  $rh_saturation->{$gene_ap}{$subtype}->[1] || 0;
               my $info =  $rh_saturation->{$gene_ap}{$subtype}->[0] || 0;
               my $cocites =  $rh_saturation->{$gene_ap}{$subtype}->[2] || 0;
                
                $sat = 1 if $cocites > 20;
                print $fh_marker_genes "$gene_ap\t$gene\t$specific_type\t$subtype\t$type\n" if $cocites > 20 ;
                $saturation{$gene_ap}{$expanded_subtype} = [$info,$sat,$cocites] ;
             }# end foreach ap cluster gene 
          } # end if saturated
       }# end gene to cluster genes
   } # end cluster
  
   $exclude{$expanded_subtype}++ if $exclude == 1;
 }# end foreach specific type
 
 return (\%saturation,\@subtypes,\%exclude,\%DIST_to_GIT);
                                                               
 
}

sub subset_ratio_network {

 my $rh_patient2subtypescore = shift;
 my $ra_patients = shift;
 my $ra_subtypes = shift;
 my $type = shift;
 my $result_dir = shift;
 my $rh_DIST_to_GIT = shift;
 my $rh_code2name = shift;
 my $rh_name2code = shift;
 my $rh_immune_cell_net = shift;
 my $case = shift;
 my $logratio_limit = shift ;
 
 my $network_dir = "$result_dir\/networks";
 unless(-e $network_dir or mkdir $network_dir) {
   die "Unable to create $network_dir\n";
 }
 
 my %patient2subtypescore = %$rh_patient2subtypescore ;
 my %mmune_cell_net = %$rh_immune_cell_net; 
 my %patient2ratio_tag = ();
 
 my %ratios = ();
 my @ratios = ();
 my %subtype_ratio2scores = ();
 my $p= 0;
 my %patient_ratio_networks = ();
 foreach my$patient(@$ra_patients){
   $p++;
   my $t = 0;
   foreach my $subtype_outer(@$ra_subtypes){
      my $ra_data = $patient2subtypescore{$patient}{$subtype_outer}; 
      my $subtype_score_outer = $ra_data->[0];
      
      foreach my $subtype_inner(@$ra_subtypes){
         next if $subtype_outer eq $subtype_inner ;
         if ($type eq 'specific') {
            $subtype_outer =~/(EXT_[a-zA-Z0-9]+_)/;
            my $general_1  = $1 ;
            $subtype_inner =~ /(EXT_[a-zA-Z0-9]+_)/;
            my $general_2  = $1 ;
            next if $general_1 eq $general_2 ;
            next if $subtype_score_outer == 0;
         }

         my $ra_data_2 = $patient2subtypescore{$patient}{$subtype_inner}; #[$data_point,$q];
         my $subtype_score_inner = $ra_data_2->[0] ; 
         
         next if $subtype_score_inner == 0;
         my $ratio =  ($subtype_score_outer / $subtype_score_inner ) ;
         next if $ratio < 0;
         my $log_ratio = log2($ratio); 
         
         my $patient_tag = $patient.'_'.$subtype_outer.'_'.$subtype_inner;
         $patient2ratio_tag{$patient_tag} = $ratio;
         
         $t++;
         if ($log_ratio >= $logratio_limit) {
            # construct the patient specific network on the basis of the edge information for DIST only.
            if ($type eq 'specific') {
               my $general_outer = $rh_DIST_to_GIT->{$subtype_outer};
               my $general_inner = $rh_DIST_to_GIT->{$subtype_inner};
               push(@{$patient_ratio_networks{$patient}}, [$general_outer,$general_inner,$subtype_outer,$subtype_inner,$log_ratio]); 
            }
            
         }
         my $ratio_tag = $subtype_outer.'_'.$subtype_inner;
         push (@ratios, $ratio_tag)  unless exists  $ratios{$subtype_outer}{$subtype_inner};
         push (@{$subtype_ratio2scores{$ratio_tag}}, [$ratio,$patient]);
         
         $ratios{$subtype_outer}{$subtype_inner}++;
         $ratios{$subtype_inner}{$subtype_outer}++; 
      }
   }
   
 }# end for each patient
 
   my $network_file_code = "$result_dir\/$case\_immunecell_net_codes.txt";
   my $fh_net_code = new IO::File "> $network_file_code" or die "cannot open $network_file_code...$!";
   my $i = 0;
   my %patient_edge_count = ();
   my %patient_edge_covered = ();
   my $pn = scalar(@$ra_patients) / 10; 
   foreach my$patient(@$ra_patients){
      
      my $ra_network_edges = $patient_ratio_networks{$patient};
      next if !defined($ra_network_edges); 
      my %edge_covered = ();
      
      # network related files
      my $network_file = "$network_dir\/$patient\_immunecell_net.txt";
      my $network_edges = "$network_dir\/$patient\_immunecell_associated_DIST_edges.csv";
   
      my $fh_net = new IO::File "> $network_file" or die "cannot open $network_file...$!";
      my $fh_net_edges = new IO::File "> $network_edges" or die "cannot open $network_edges...$!";
      my $e=0;
      foreach my $ra_edge(@$ra_network_edges){
         my $git1 = $ra_edge->[0];
         my $git2 = $ra_edge->[1];
         my $specific_type1 = $ra_edge->[2];
         my $specific_type2 = $ra_edge->[3];
         my $edge_score = $ra_edge->[4]; # log2 ratio
         
         print $fh_net_edges "$patient,$git1,$git2,$specific_type1,$specific_type2,$edge_score\n";
         my $code1 = $rh_name2code->{$git1};
         my $code2 = $rh_name2code->{$git2}; 
         
         my $i =  $patient_edge_count{$git1}{$git2};
         unless (exists $edge_covered{$git1}{$git2}){
            $e++;
            print $fh_net "$git1\t$git2\tpd\t$edge_score\n"; 
            if (exists $mmune_cell_net{$code1}{$code2} && $i > $pn) {
                print $fh_net_code "$code1\t$code2\tpd\t$i\n" unless (exists $patient_edge_covered{$code1}{$code2});
                $patient_edge_covered{$code1}{$code2}++;
                $patient_edge_covered{$code2}{$code1}++;
            } 
            $patient_edge_count{$git1}{$git2}++;
            $patient_edge_count{$git2}{$git1}++;
         }
         $edge_covered{$git1}{$git2}++;
         $edge_covered{$git2}{$git1}++;
      }
   } 
 return (\%patient2ratio_tag, \@ratios,\%subtype_ratio2scores);
}


sub cell_specific_networks{
   
my %cell_specific_networks = (
                              'cytotoxicT' => [
                                               'NaiveCD8+T-cell_Ce',
                                               'CD8+CentralMemory_Ce',
                                               'CD8+EffectorMemory_Ce',
                                               'EffectorMemoryRA',
                                               'HumanChronicHIVSampleCD8+Tcel',
                                               'HumanAcuteHIVSampleCD8+Tcel',
                                               'HumanNon-progressorHIVSampleCD8+Tcel',
                                               'HumanUninfectedSampleCD8+Tcel',
                                               'Tet+Tre',
                                               'Tet-Tre',
                                               'NaiveCD8+CD3+Tcel',
                                               'PD-1highCD8+CD3+Tcel',
                                               'PD-1lowCD8+CD3+Tcel',
                                               'IE-CTLNKG',
                                               'PB-CTLNKG2',
                                               'HumanChronicHIVCD4+Tcel',
                                               'HumanNon-progressorHIVCD4+Tcel',
                                               'HumanAcuteHIVCD4+Tcel',
                                               'HumanChronicHIVCD8+Tcel',
                                               'HumanAcuteHIVCD8+Tcel',
                                               'HumanNon-progressorHIVCD8+Tcel',
                                               'HumanUninfectedCD8+Tcel',
                                               'CD8+T-cellsfromperipheralblo',
                                               'CD8+TcellsmR',
                                               'BloodCD8Temcel',
                                               'IntraepithelialCD8Temcel',
                                               'CD8_T_cell_pool_melano',
                                               'CD8_T_cell_pool_healt'
                                               ],
                                               
                              'natural_killer' => [
                                                   'MatureNKcell_CD56-CD16-CD3-_Ce',
                                                   'MatureNKcell_CD56+CD16+CD3-_Ce',
                                                   'MatureNKcell_CD56-CD16+CD3-_Ce',
                                                   'NKcells_unstimulat',
                                                   'Resting_CD8Tcel',
                                                   'NK_time2hou',
                                                   'Resting_',
                                                   'NK_timepoint8hou',
                                                   'NK_timepoint24hou',
                                                   'Lymphoidcontr',
                                                   'NK_Donor',
                                                   'ENK_myeloma_patient',
                                                   'ENK_Donor',
                                                   'NK_Donor',
                                                   'CD56+NK-cellsfromperipheralblo',
                                                   'NKcellsmR',
                                                   'NK_cell_pool_healt',
                                                   'NK_cell_pool_melano'
                                                  ],
                                                   
                              'Th1' => [
                                        ' NaiveCD4+T-cell_Ce',
                                        'CD4+EffectorMemory_Th1',
                                        'CD4+CentralMemory_Th1',
                                        'NaiveCD4+T-cell_Th1',
                                        'NaiveCD4+T-cell_Ce',
                                        'HumanChronicHIVSampleCD4+Tcel',
                                        'HumanAcuteHIVSampleCD4+Tcel',
                                        'HumanNon-progressorHIVSampleCD4+Tcel',
                                        'HumanUninfectedSampleCD4+Tcel',
                                        'EffectormemoryTcells_unstimulat',
                                        'CentralmemoryTcells_unstimulat',
                                        'Th1cells_unstimulat',
                                        'CD4+T-cellsfromperipheralblo',
                                        'TIL', # 0 hours 
                                        'TIL2',# 24 hours 
                                        'RO',# 0 hours ?
                                        'RO2',# 24 hours
                                        'PBDonor',
                                        'TILPatient',
                                        'LNPatient',
                                        'PBPatient' ,
                                        'CD4+TcellsmR',
                                        'BloodCD4Temcel',
                                        'LaminapropriaCD4Temcel',
                                        'IntraepithelialCD4Temcel',
                                        'CD4_T_cell_pool_melano',
                                        'CD4_T_cell_pool_healt'
                                       ],
                                     
                              'Th2' => [
                                        ' NaiveCD4+T-cell_Ce',
                                        'CD4+EffectorMemory_Th2',
                                        'CD4+CentralMemory_Th2',
                                        'NaiveCD4+T-cell_Th2',
                                        'HumanChronicHIVSampleCD4+Tcel',
                                        'HumanAcuteHIVSampleCD4+Tcel',
                                        'HumanNon-progressorHIVSampleCD4+Tcel',
                                        'HumanUninfectedSampleCD4+Tcel',
                                        'EffectormemoryTcells_unstimulat',
                                        'CentralmemoryTcells_unstimulat',
                                        'Th2cells_unstimulat',
                                        'CD4+T-cellsfromperipheralblo',
                                        'TIL', # 0 hours 
                                        'TIL2',# 24 hours 
                                        'RO',# 0 hours ?
                                        'RO2',# 24 hours 
                                        'PBDonor',
                                        'TILPatient',
                                        'LNPatient',
                                        'PBPatient',
                                        'CD4+TcellsmR',
                                        'BloodCD4Temcel',
                                        'LaminapropriaCD4Temcel',
                                        'IntraepithelialCD4Temcel',
                                        'CD4_T_cell_pool_melano',
                                        'CD4_T_cell_pool_healt'
                                       ],
                                       
                              'B_cells' => [
                                           'EarlyB-cell_Ce',
                                           'MatureB-cellclassabletoswitch_Ce',
                                           'MatureB-cellclassswitched_Ce',
                                           'MatureB-cells_Ce',
                                           'NaiveB-cell',
                                           'ProB-cell_Ce',
                                           'Bcells_unstimulat' ,
                                           'CD19+B-cellsfromperipheralblo',
                                           'BcellsmR',
                                           'B_cell_pool_melano',
                                           'B_cell_pool_health',
                                           'NaïveB-cells_Ce'
                                          ],
                              
                              'Dendritic' => [
                                              'MyeloidDendriticCell',
                                              'Immaturedendriticcel',
                                              'Dendriticcells_LPS4',
                                              'HumanDendriticCellL.major24hrsbiologicalre',
                                              'HumanDendriticCellL.major8hrsbiologicalre',
                                              'HumanDendriticCellL.major4hrsbiologicalre',
                                              'HumanDendriticCellL.major2hrsbiologicalre',
                                              'HumanDendriticCelluninfected2hrsbiologicalre' ,
                                              'mDCmR',
                                              'pDCmR',
                                              'MyeloidDendriticCell_Ce',
                                              'PlasmacytoidDendriticCell_Ce'
                                             ],
                              
                              'treg' => [
                                         'Treg',
                                         'TH' ,
                                         'PCa1Tr',
                                         'PCa2Tr',
                                         'PCa3Tr',
                                         'HD1Tr',
                                         'HD2Tr',
                                         'HD3Tr',
                                         'F-S-',
                                         'HI_F+S+ctl_freshA',
                                         'EFF_F-S-ctl_fre',
                                         'FEEDERCTL2A',
                                         'F+S+',
                                         'F+S-ACT',
                                         'F+S+ACT',
                                         'FEEDERCTL1A',
                                         'FEEDERCTL4A',
                                         'EFF_F-S-ctl_freshA',
                                         'HI_F+S+ctl_fre',
                                         'F+S-',
                                         'FEEDERCT',
                                         'F-S-ACT'
                                        ],
                                         
                              'Macrophage' => [ 
                                               'Monocyte',
                                               'Monocyteat',
                                               'classicalorM1activatedmacrophag',
                                               'AlternativeorM2activatedmacrophag',
                                               'Monocyteat3da',
                                               'Macrophageat7da',
                                               'Macrophages_unstimulat',
                                               'Macrophages_LPS',
                                               'Primaryhumanmonocytesfromhealthycontrol',
                                               'Primaryhumanmonocytesfrompatientwithgram-negativesepsis',
                                               'Primaryhumanmonocytesfrompatientwithmetastaticbreastcancer',
                                               'Primaryhumanmonocytesfrompatientwithtuberculosis',
                                               'DcR3-treatedMDMs_2days_re',
                                               'hIgG1-treatedMDMs_2days_re',
                                               'PBMo-derivedmacrophag',
                                               'PBMo-derivedmacrophagesstimulatedwithTE-9',
                                               'PBMo-derivedmacrophagesstimulatedwithTE-8',
                                               'PBMo-derivedmacrophagesstimulatedwithTE-15',
                                               'MY',
                                               'MEP',
                                               'BC',
                                               'early_PM',
                                               'late_PM',
                                               'GMP',
                                               'MM',
                                               'PMN',
                                               'HSC',
                                               'MPP',
                                               'CMP',
                                               'CD14+monocytesfromperipheralblo',
                                               'HD_MONOCYT',
                                               'MONOCYT',
                                               'MonocytesmR',
                                               'Commonmyeloidprogenitor_Ce',
                                               'ColonyFormingUnit-Monocyte_Ce',
                                               'Monocyte_Ce'
                                               ],
                                               
                             'Granulocytes' => [
                                                'Eosinophils_PMA_',
                                                'Eosinophils_contr',
                                                'Basophils_unstimulat',
                                                'Neutrophils_unstimulat',
                                                'Neutrophils_LPS',
                                                'Basophils_unstimulat',
                                                'Neutrophils_unstimulat',
                                                'Neutrophils_LPS',
                                                'MSC_B_',
                                                'Resto_C_',
                                                'Resto_C_wo',
                                                'MSC_A_',
                                                'MSC_B_wo',
                                                'Resto_B_wo',
                                                'Resto_B_',
                                                'CD15+neutrophilsfromperipheralblo',
                                                'NeutrophilsmR',
                                                'EosinophilsmR',
                                                'Granulocyte_NeutrophilicMetamyelocyte_Ce',
                                                'Granulocyte_monocyteprogenitor_Ce',
                                                'ColonyFormingUnit-Granulocyte_Ce',
                                                'Granulocyte_Neutrophil__Ce',
                                                'Eosinophill_Ce',
                                                'Basophils_Ce'
                                               ],
                            'Mast_cells' => [
                                             'Cordblood-derivedmastcells_IgE',
                                             'Cord_blood_derivedmast_cells_contr'
                                             ]
                             
                              );

 return \%cell_specific_networks;
}


sub fetch_saturation {

 my $ra_subtypes = shift;
 my $folder = shift;
 my $fh_marker_genes_GIT = shift;
 my $saturation_limit = shift;

 my %gene2subtype = ();
 my %subtype2genes = ();
 my $count = 0;
  foreach my$subtype(@$ra_subtypes){
     my $file = "$folder\/$subtype\_res.csv";
     #print "Creating file handle for $file\n";
     my $fh = new IO::File $file or die "cannot open $file...$!";
     my $c=0;
     while (! $fh->eof()){
        my $line;
        while (defined ($line = <$fh>)){
           $c++;
           $line  =~s/\s+$//;
           my @row = split("\,", $line);
           my $sym = $row[1];  $sym  =~s/\s+$//;
           my $information = $row[2];  next if (!defined($information));
           $information  =~s/\s+$//;
           my $saturation = $row[3];  next if (!defined($saturation));
           my $cocites = $row[4];  next if (!defined($cocites));
           $saturation  =~s/\s+$//;
           
           if ($subtype eq 'Granulocytes') {
            next if $sym eq 'CD8';
            next if $sym eq 'CD8A';
            next if $sym eq 'CD4';
            next if $sym eq 'IFNG';
            next if $sym eq 'TBET';
            next if $sym eq 'IL4';
            next if $sym eq 'GATA3';
            next if $sym eq 'FOXP3';
            next if $sym eq 'IL2RA';
            next if $sym eq 'TGFB1'; 
            next if $sym eq 'IL5'; 
            next if $sym eq 'IL3';
            next if $sym eq 'IL2'; 
            next if $sym eq 'IL13'; 
            next if $sym eq 'TBX21'; 
            next if $sym eq 'TNF'; 
            next if $sym eq 'CSF2';
            next if $sym eq 'IL8'; 
            next if $sym eq 'CD83';
            next if $sym eq 'ITGAX';
            next if $sym eq 'CD1A'; 
            next if $sym eq 'CCL22'; 
            next if $sym eq 'IL17'; 
            next if $sym eq 'IL17A'; 
            next if $sym eq 'IL17B'; 
            next if $sym eq 'IL17C'; 
            next if $sym eq 'IL17D'; 
            next if $sym eq 'IL17E'; 
            next if $sym eq 'IL17RA';
            next if $sym eq 'IL17RB'; 
            next if $sym eq 'IL17RC'; 
            next if $sym eq 'IL17RD'; 
            next if $sym eq 'IL17RE'; 
            next if $sym eq 'IL22'; 
            next if $sym eq 'CD56';  
            next if $sym eq 'NCAM'; 
            next if $sym eq 'CD16'; 
            next if $sym eq 'MS4A1';
            next if $sym eq 'CD19'; 
            next if $sym eq 'ITGAM';

            next if $sym eq 'INS';
            next if $sym eq 'FAM48A'; 
            next if $sym eq 'JUN' ; 
            next if $sym eq 'MS' ; 
            next if $sym =~/MAPK/; 
            next if $sym eq 'CD80'; 
            next if $sym eq 'CD86'; 
            next if $sym eq 'CD40';
            next if $sym eq 'CD40LG';
            next if $sym eq 'CD1D';
            next if $sym eq 'CD14';
            next if $sym eq 'CCR7';
            next if $sym eq 'NDUFA2';
            next if $sym eq 'ITGB2';
            next if $sym eq 'FASLG';
    
            next if $sym eq 'IL18'; 
            next if $sym eq 'IFNA1'; 
            next if $sym eq 'IL15';  
            next if $sym eq 'IL2';  
            next if $sym eq 'IL10';  

           } 
           
           $count++ if $saturation >=  $saturation_limit;
           print $fh_marker_genes_GIT "$count\t$sym\t$subtype\n" if $saturation >= $saturation_limit;
           $gene2subtype{$sym}{$subtype} = [$information,$saturation,$cocites] ;
           push(@{$subtype2genes{$subtype}}, $sym); 
        }
      } 
     
  }
 return \%gene2subtype,\%subtype2genes;
}


sub log2 {
   my $n = shift;
   if ($n == 0) {
      return 0;
   }
   
   return (log($n)/log(2));
}




1;
