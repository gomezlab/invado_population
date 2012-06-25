#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use Image::ExifTool;
use Getopt::Long;
use Data::Dumper;

use Config::Adhesions;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "queue=s", "resource|R=s") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################

my $exp_type_folder = "'".dirname($cfg{exp_results_folder}) ."'";

my @matlab_code;
push @matlab_code, "build_raw_montage_file_set($exp_type_folder,'gel.png')";
push @matlab_code, "build_raw_montage_file_set($exp_type_folder,'puncta.png')";
push @matlab_code, "build_raw_montage_file_set($exp_type_folder,'uncorrected_puncta.png')";
push @matlab_code, "build_raw_montage_file_set($exp_type_folder,'uncorrected_gel.png')";
push @matlab_code, "build_raw_montage_file_set($exp_type_folder,'highlighted_mask.png')";
push @matlab_code, "build_raw_montage_file_set($exp_type_folder,'cell_overlaps.png')";
push @matlab_code, "build_raw_montage_file_set($exp_type_folder,'invader_vis.png')";
push @matlab_code, "build_raw_montage_file_set($exp_type_folder,'gel_highlights.png')";
push @matlab_code, "build_raw_montage_file_set($exp_type_folder,'tracking.png')";

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'montage_vis');
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################
