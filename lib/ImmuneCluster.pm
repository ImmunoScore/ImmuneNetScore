##############################################################################
# file: ImmuneCluster.pm
# desc: Clustering usng Affinity propagation
##############################################################################
package ImmuneCluster;

use IO::File;
use strict;


sub specific_type_apclusters {

 my $DIR = shift;
 
 opendir(Dir, $DIR) or die "cannot open directory $DIR $!";
 my @files = grep(/\.txt$/,readdir(Dir));	
 
 my $f = 0;

 my %specific_type2gene_apclusters= ();
 my %specific_type2gene_cluster_IDs= ();
 foreach my$file(@files){
   $f++;
   $file  =~/apclusters_(.+)\.txt/;
   my $specific_type = $1;
   my $fh = new IO::File "$DIR\/$file" or die "cannit oen file $DIR\/$file $!";
   
   my @cluster_IDs = ();
   while (! $fh->eof()){
   my $line;
   while (defined ($line = <$fh>)){
      $line  =~s/\s+$//;
      my @line = split("\t",$line);
      my $no = $line[0];
      my $cluster_ID = $line[1];
      my $cluster_genes = $line[2];
      my @genes = split(",", $cluster_genes);
      push (@{$specific_type2gene_cluster_IDs{$specific_type}},$cluster_ID); 
      my %gene2clustergenes = ();
      foreach my$gene_outer(@genes){
	 foreach my$gene_inner(@genes){
	    next if $gene_inner eq $gene_outer;
	    push (@{$gene2clustergenes{$gene_outer}},$gene_inner ); 
	 }
      }
      $specific_type  =~s/\s+$//; 
      $specific_type2gene_apclusters{$specific_type}{$cluster_ID} = \%gene2clustergenes; 
      }  
   }

 }
   
 return (\%specific_type2gene_apclusters,\%specific_type2gene_cluster_IDs); 
}


1;
