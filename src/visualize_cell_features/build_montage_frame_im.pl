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

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{convert_to_bmp} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "i_name=s", "convert_to_bmp") or die;

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
my @field_order = (
	4,5,14,15,24,
	3,6,13,16,23,
	2,7,12,17,22,
	1,8,11,18,21,
	0,9,10,19,20,
	);

@fields = @fields[@field_order];

my @image_dirs = <$fields[0]/individual*/*>;
my $image_count = scalar(@image_dirs);

foreach my $i_num (0..($image_count - 1)) {
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
	
	my $extra;

	my $command = "montage " . join(" ", @file_list) . " -geometry x160+0+0 -tile 5x5 $out_file_png";
	if ($opt{debug}) {
		print "$command\n";
	} else {
		system($command);
	}

	if (@level_norms) {
		my $command = "convert $out_file_png -level $level_norms[0],$level_norms[1] $out_file_png;";
		if ($opt{debug}) {
			print "$command\n";
		} else {
			system($command);
		}
	}
	
	#Add field number annotation to the first image
	if ($i_num == 0) {
		my $annotate_set = "";
		my $x_gap = 210;
		my $y_gap = 160;
		my @x_pos = (5 + $x_gap*0,5 + $x_gap*1,5 + $x_gap*2,5 + $x_gap*3,5 + $x_gap*4) x 5;
		my @y_pos = ((15+$y_gap*0) x 5, (15+$y_gap*1) x 5, (15+$y_gap*2) x 5,
			(15+$y_gap*3) x 5, (15+$y_gap*4) x 5);
		for my $i (0..$#field_order) {
			my $field_num = $field_order[$i] + 1;
			$annotate_set .= "-annotate +$x_pos[$i]+$y_pos[$i] $field_num ";
		}
		
		my $command = "convert $out_file_png -pointsize 16 $annotate_set $out_file_png;";
		if ($opt{debug}) {
			print "$command\n";
		} else {
			system($command);
		}
	}
	
	my $out_file_bmp = "$output_folder/$opt{i_name}/" . sprintf("%02d",$i_num+1) . ".bmp";
	if ($opt{convert_to_bmp}) {
		my $command = "convert $out_file_png $out_file_bmp; rm $out_file_png;";
		if ($opt{debug}) {
			print "$command\n";
		} else {
			system($command);
		}
	}
}
