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
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################

my @image_folders = sort <$cfg{individual_results_folder}/*>;
die "Expected to find results directories." if (scalar(@image_folders) == 0);

if ($opt{debug}) {
    if (scalar(@image_folders) > 1) {
        print "Individual results folders found: $image_folders[0] - $image_folders[$#image_folders]\n";
    } elsif (scalar(@image_folders) == 0) {
        warn "Couldn't find any results folders in $cfg{individual_results_folder} subfolders\n\n";
    } else {
        print "Results folder found: $image_folders[0]\n";
    }
}

my @matlab_code = &create_all_matlab_commands(@image_folders);

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'cell_mask_props');
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub create_all_matlab_commands {
	my @image_folders = @_;
    my @matlab_code;

    foreach (0..$#image_folders) {
		if ($_ == 0) {
			$matlab_code[0] .= "find_cell_mask_properties('$image_folders[0]','$image_folders[0]')\n";
		} else {
			$matlab_code[0] .= "find_cell_mask_properties('$image_folders[$_]','$image_folders[$_ - 1]')\n";
		}
    }

    return @matlab_code;
}
