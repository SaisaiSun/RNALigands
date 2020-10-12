#!/usr/bin/perl -w
use strict;

###################################
#This code is used to run RNAstructure using Fold algrithom to predict Secondary structure for RNA.

#Saisai Sun
#Date:6-5-2020

##################################


my $query_id=$ARGV[0];
my $input_dir=$ARGV[1];
my $output_dir=$ARGV[2];
#my $RNAstr_dir="/home/sssun/software/RNAstructure/exe";
my $RNAfold_dir="/var/www/rnaligands/ViennaRNA/bin";

#`$RNAstr_dir/Fold $input_dir/$query_id.fasta $output_dir/$query_id.ct --loop 30 --maximum 20 --percent 10 --temperature 310.15 --window 3`;
#`$RNAstr_dir/ct2dot $output_dir/$query_id.ct 1 $output_dir/$query_id\_dot.txt`;
`$RNAfold_dir/RNAfold --noPS -i $input_dir/$query_id.fasta >$output_dir/$query_id\_dot.txt`;

