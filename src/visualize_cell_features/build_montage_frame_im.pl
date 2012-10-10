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
use Getopt::Long;
use Data::Dumper;
use Text::CSV;

use Config::Adhesions;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "i_name=s") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};
die "Can't find file name to montage, expected i_name 'filename'" if not exists $opt{i_name};

#strip off the .png if specified on the i_name parameter
$opt{i_name} =~ s/\.png$//;

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################

#determine level norms, if they exist
my $norm_file;
if ($opt{i_name} eq "gel") {
	$norm_file = "$cfg{exp_results_folder}/../all_field_resources/gel_image_range.csv";
} elsif ($opt{i_name} =~ /puncta/) {
	$norm_file = "$cfg{exp_results_folder}/../all_field_resources/puncta_image_range.csv";
}

my @level_norms;
if (defined $norm_file) {
	open CSV, "$norm_file" or die "$!";
	my @norm_lines = <CSV>;
	close CSV;

	my $csv = Text::CSV->new();
	
	$csv->parse($norm_lines[1]);
	@level_norms = $csv->fields();
}

#build/execute the montage image commands
my $output_folder = "$cfg{exp_results_folder}/../montage/";
mkpath($output_folder);

my @fields = <$cfg{exp_results_folder}/../time_series*>;
@fields = @fields[
	4,5,14,15,24,
	3,6,13,16,23,
	2,7,12,17,22,
	1,8,11,18,21,
	0,9,10,19,20,
	];

my @image_dirs = <$fields[0]/individual*/*>;
my $image_count = scalar(@image_dirs);

foreach my $i_num (0..($image_count - 1)) {
# foreach my $i_num (0) {
	my @file_list;
	foreach my $this_field (@fields) {
		my @images = <$this_field/individual*/*/$opt{i_name}.png>;
		if (! -e $images[$i_num]) {
			print "Can't find $images[$i_num] $this_field\n";
		} 
		push @file_list, $images[$i_num];
	}
	
	mkpath("$output_folder/$opt{i_name}");
	
	my $out_file_png = "$output_folder/$opt{i_name}/" . sprintf("%02d",$i_num+1) . ".png";
	my $out_file_bmp = "$output_folder/$opt{i_name}/" . sprintf("%02d",$i_num+1) . ".bmp";
	
	my $extra;

	my $command = "montage " . join(" ", @file_list) . " -geometry x160+0+0 -tile 5x5 $out_file_png";
	
	if ($opt{debug}) {
		print "$command\n";
	} else {
		system($command);
	}
	
	if (@level_norms) {
		my $command = "convert $out_file_png -level $level_norms[0],$level_norms[1] $out_file_png; convert $out_file_png $out_file_bmp;";
		if ($opt{debug}) {
			print "$command\n";
		} else {
			system($command);
		}
	} else {
		my $command = "convert $out_file_png $out_file_bmp";
		if ($opt{debug}) {
			print "$command\n";
		} else {
			system($command);
		}
	}
	my $command = "rm $out_file_png";
	if ($opt{debug}) {
		print "$command\n";
	} else {
		system($command);
	}
}
