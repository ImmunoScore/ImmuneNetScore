#############################################################################
# file: ImmuneGEO.pm
# desc: functions which interlink with the GEO data in ImmuneScore
##############################################################################

package ImmuneGEO;

use strict;
use IO::File;


sub fetch_platform_annotations {

 my $platform_DIR = shift;
 my $name = shift;
 my $begin = shift;
 my $sym_index = shift;
 
 my %illumina2gene= (); 
 
 # source GEO reference 
 my $fh = new IO::File "$platform_DIR\/$name";
 
 my $i = 0;
 my $sym_delim='///';
 while (! $fh->eof()){
   my $line;
   while (defined ($line = <$fh>)){
     $line  =~s/\s+$//; $line  =~s/"//g;
     $i++;
     next if $i < $begin;
     my @row = split("\t", $line);
     my $id = $row[0];  $id  =~s/\s+$//;
     my $syms = $row[$sym_index] || "-";
     
     if($syms =~/\s+(\/\/)\s+/){
      $sym_delim = $1;
      print "$sym_delim\n"; 
     }
     
     my @syms = split("$sym_delim",$syms );
        
     my $s = 0;
     foreach my$sym(@syms){
	  $s++;
	  next if (exists $illumina2gene{$id});
	  
	  $sym  =~s/\s+$//;
	  $sym  =~s/^\s+//;  
	  $illumina2gene{$id} = $sym if ($s == 1 && $sym_delim eq '///');
	  $illumina2gene{$id} = $sym if ($s == 2 && $sym_delim eq '//');
     }
     
   }
 }
 
 return \%illumina2gene;
 
}


sub config_immune {

 my $fh = shift;
 
 my %config= ();
 my %properties = ();

 my $property_type='';
 my $values='';
 my $regexp='';
 my $patient_color = '';
 while (! $fh->eof()){
   my $line;
   while (defined ($line = <$fh>)){
     $line  =~s/\s+$//; $line  =~s/"//g;

     if($line  =~/^#{1}(\w+):\s*(.+)/){
	my $arg = $1;  $arg =~s/\s+$//;$arg =~s/^\s+// ;
	my $arg_value = $2; $arg_value =~s/\s+$//; $arg_value =~s/^\s+// ;
        #print "$arg\t$arg_value\n";
	$config{$arg} = $arg_value
     }
     
    # for config file with patient groups and clinical properties
    if($line  =~/^\/{1}(property_type):\s*(.*)/){
        $property_type = $2; $property_type =~s/\s+$//; $property_type =~s/^\s+// ;
     }
    
    if($line  =~/^\/{1}(properties):\s*(.*)/){
        $values = $2;  $values =~s/\s+$//; $values =~s/^\s+// ;
     }
    
    if($line  =~/^\/{1}(regexp):\s*(.*)/){
        $regexp = $2 ;
     }
    
    if($line  =~/^\/{1}(patient_row_color):\s*(.*)/){
        $patient_color = $2 ;
    }
    
    if($line  =~/^\/{2}/){
	#print "$property_type\t$values\t$regexp\n";
	$property_type =~s/\s+$//; $property_type =~s/^\s+// ;
	$property_type  =~/(.+)\=(.+)/;  
        my $assigned_property = $1;
        my $GEO_annatation = $2;
	$GEO_annatation =~s/\s+$//; $GEO_annatation =~s/^\s+// ;
	$assigned_property =~s/\s+$//; $assigned_property =~s/^\s+// ;
	#print "$GEO_annatation, $assigned_property\t$values\t$regexp\n";
	my $ra_entry = parse_properties($values,
					$GEO_annatation,
					$regexp,
					$patient_color);
	
	$properties{$assigned_property} = $ra_entry;
	$values='';
	$property_type='';	
     }# 
      
   }
 }

 return(\%config,\%properties); 
 
}

sub parse_properties{
 my $values = shift ;
 my $GEO_annatation = shift;
 my $regexp = shift;
 my $patient_color = shift;
 
 my @values = split (",", $values);
 my @assigend_properties=();
 my %annotated_properties = ();
 my %property = ();
 foreach my $value(@values){
    $value  =~/(.+)\=(.+)/; $value =~s/\s+$//; $value =~s/^\s+// ;
    my $assigned_value = $1; $assigned_value =~s/\s+$//; $assigned_value =~s/^\s+// ;
    my $GEO_value = $2; $GEO_value =~s/\s+$//; $GEO_value =~s/^\s+// ;
    
    push(@assigend_properties,$assigned_value);
    $annotated_properties{$assigned_value} = $GEO_value;
  
 }
 
 return [\%property,
	 \@assigend_properties,
	 $GEO_annatation,
	 \%annotated_properties,
	 $regexp,
	 $patient_color 
	 ]; 
 
}

1;
