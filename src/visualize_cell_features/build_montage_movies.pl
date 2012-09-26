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

my $montage_folder = catdir(dirname($cfg{exp_results_folder}), 'montage');

# my $base_movie_command = "ffmpeg -sameq -v 0 -y -r $cfg{movie_frame_rate} ";
my $base_movie_command = "avconv -y -r $cfg{movie_frame_rate} ";

my @movie_folders = <$montage_folder/*>;
@movie_folders = grep -d $_, @movie_folders;

for my $this_montage_folder (@movie_folders) {
	my @files = <$this_montage_folder/*>;
	my $image_num_length = length(@files);

	my $movie_output_file = catfile($montage_folder,'../',basename($this_montage_folder).".mp4");
	
	my $movie_command = $base_movie_command . 
		"-i $this_montage_folder/%0" . $image_num_length . "d.png -q 1 $movie_output_file ";
	
	if ($opt{debug}) {
		print $movie_command, "\n";
	} else {
		system("$movie_command");
	}
}

