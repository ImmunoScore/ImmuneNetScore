# ImmunoScore
A computational framework to profile networks of distinct immune-cells from tumor transcriptomes

ImmuneScore VERSION 1.0            User  Perl Documentation           
NAME
       ImmunoScore.pl  

SYNOPSIS
       Reads in an expression profile for a quantile normalized datasets in the NCBI
       GEO data base format and calculates an immune profile for both general
       immune-cell types (GITs) and distinct immune-cell types (DISTs) based on
       an integrated framework of knowledge based information for human genes 

       Usage:
       
       $perl ImmunoScore.pl [CASE_NAME] [CONFIG_FILE] [GEO_SERIES_MATRIX_FILE] 

DESCRIPTION

       The script analyzes a transcriptome of a complex biological sample, such as a tumor,
       and calculates immune profiles for both general immune-cell types (GITs) and distinct
       immune-cell types (DISTs) based on an integrated framework of knowledge based information
       for human genes. More detailed background information is accessible in the manuscript:
       A computational framework to profile networks of distinct immune-cells from tumor transcriptomes,
       Clancy and Hovig, 2015. 
       
       The script accesses key sourced data from './data/' directories in the repository neccessary to calculate 
       immune-cell profiles and networks of immune-cell  whose signatures expression profiles have been inferred from
       transcriptome analysis. The following are the source data in subdirectories  './data/': 
          - saturation: GIT profiles with gene saturation scores computed from all of Medline
          - aplcusters: The resltant affinity propagation clusters computed from DIST specific
                        protein-protien interaction networks
          - platform: Expression array plaform data from the GEO_series_matrix_file
       
       The script then creates a results directory with the following .TXT files detailing the GIT
       and DIST immmune profiles for the transcritpome being analyzed in './CASE_NAME':
          - general_res_CASE_NAME_subtype_scores.txt - calculated GIT scores for all samples/patients
          - specific_res_CASE_NAME_subtype_scores.txt - calculated DIST scores for all samples/patients
       
       The script creates a directory for patient/sample specific immune-cell networks in
       the './CASE_NAME/networks' directory containing the following file types:
          - GSMXXXX_immunecell_net.txt - network files
          - GSMXXXX_immunecell_associated_DIST_edges.csv - underlying DIST connections for the network
      

INPUT
       -CASE_NAME
          Assigned name to the study/case. A '../CASE_NAME/' directory is created, score files and
          '../CASE_NAME/results' sub-directories
          
       - CONFIG_FILE [OPTIONS]
          Each option in the configuration files has a prefix '#' followed by the option and its variabe: '#OPTION: VAR'
          - #begin
          - #case: CASE_NAME
          - #accession: GSEnnnn
          - #platform: GPLmmmm
          - #gene_platform:  Row number from GPLmmmm indicating start of ID row in the platform
          - #begin_data_platform: Row number from GPLmmmm indicating start of data row in the platform
          - #patient_id_row: Row number from GSEnnnn indicating start of patient row in the acession
          - #sample_id_row:  Row number from GSEnnnn indicating start of sample ID row in the acession
          - #saturation: value between 0.0 and 1.0. The limit of literature saturation of genes to immune-cells (default=0.9)
          - #logratio: :value between 0.0 and 1.0. The log ratio which defines network connection between two immune-cells (default=0.5)
          - ##end config

       -GEO_SERIES_MATRIX_FILE
          This format is available through the NCBI Gene Expression Omnibus (GEO) website.
          The files are typically named GSEnnnn-GPLmmm_series_matrix.txt (where nnnn is the series number and mmm is the platform number.)
          These files only include samples for a single platform.
          Files downloaded from GEO are usually compressed, with .gz extensions. They need to be decompressed (e.g. using gunzip) before use. 

OUTPUT
      
       - General Immune-cell Type (GIT) scores in the filename:  general_res_CASE_NAME_subtype_scores.txt
       (calculated GIT scores for all samples/patients)
          Location
             '.data/CASE_NAME/results'
          Format
             Column 1 : GEO accession numbers GSMnnnnnnn for patients/samples
             Subsequent columns are the GIT scores 
          
       - Distinct Immune-cell SubType (DIST) scores in the filename:  specific_res_CASE_NAME_subtype_scores.txt
       (calculated GIT scores for all samples/patients)
          Location
             '.data/CASE_NAME/results'
          Format
             Column 1 : GEO accession numbers GSMnnnnnnn for patients/samples
             Subsequent columns are the DIST scores
             
       - Immune cell network files in the filenmae generated for each patient/sample: GSMXXXX_immunecell_net.txt 
           Location
             '.data/CASE_NAME/results/networks'
                 Format
             Column 1 : GIT node 1
             Column 2 : GIT node 2
             Column 3 : Edge directionality
                pd = directed
                pp = undirected
             Column 4 : Log ratio score
                > #logratio in CONFIG_FILE 

AUTHOR
        Trevor Clancy - trevor.clancy@rr-research.no

perl v5.8.8                       2015-12-12                        ImmuneScore VERSION 1.0
