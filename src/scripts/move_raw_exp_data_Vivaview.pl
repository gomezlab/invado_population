#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Spec::Functions;
use File::Basename;
use File::Path;
use File::Copy;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use Config::General qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "src=s", "target=s", "debug|d", "default_config|cfg|default_cfg=s", 
	"sample_names|sn=s") or die;
die "Can't find src folder on command line." if (not $opt{src});
die "Can't find target folder on command line." if (not $opt{target});
die "Can't find default_config on command line." if (not $opt{default_config});
die "Can't find sample_names on command line. Expected comma seperated list." if (not $opt{sample_names});

my $default_config_file_location = $opt{default_config};

my $relative_config_file;
if ($default_config_file_location =~ m#\.\./\.\./data/(.*)#) {
	$relative_config_file = $1;
} else {
	die "Unable to find ../../data/ in the default config file location: $default_config_file_location";
}

my %field_translation = (
	'GFP' => 'puncta',
	'RFP' => 'gel',
);

my %sample_translation; 

my $count = 1;
foreach (split(",", $opt{sample_names})) {
	$sample_translation{$count} = $_;
	$count++;
}

###############################################################################
# Main Program
###############################################################################

my @src_files = <$opt{src}/*>;
die "Couldn't find any files in $opt{src}." if scalar(@src_files) == 0;

#bundle each of the src files into the file_sets hash based on the position code
my %sample_sets;
my %log_files;
for (@src_files) {
	if ($_ =~ /Sample(\d+).*nd/) {
		# $1 - the sample number
		$log_files{$1} = $_;
		next;
	}

	#skip over thumbnail files
	next if ($_ =~ /thumb/); 

	if ($_ =~ /Sample(\d+).*_w\d+(.*)_s(\d+)_t(\d+).TIF/) {
		# $1 - the sample number
		# $2 - the type of image, puncta or gel, based on translations above
		# $3 - the stage position in that sample
		# $4 - the time point number
		$sample_sets{$1}{$3}{$2}[$4] = $_;
	} else {
		die "Unrecoginable filename ($_), can't find position code.";
	}
}

foreach (keys %sample_sets) {
	&confirm_stage_hash(%{$sample_sets{$_}});
}

my $data_folder = &determine_data_folder_name($opt{default_config});
foreach my $sample_num (sort keys %sample_sets) {
	my %stage_position_sets = %{$sample_sets{$sample_num}};
	my $sample_folder = $sample_translation{$sample_num};
	die "Couldn't find sample name for sample: $sample_num" if (not defined $sample_folder);
	
	foreach my $position (sort {$a<=>$b} keys %stage_position_sets) {
		my $position_directory = catdir($opt{target}, $sample_folder, 
			"time_series_" . sprintf('%02d',$position));
		
		my $exp_name;
		if ($position_directory =~ /$data_folder(.*)/) {
			$exp_name = $1;
		} else {
			die "Unable to find exp name in $position_directory";
		}

		my $log_file = $log_files{$sample_num};
		
		if ($opt{debug}) {
			print "mkpath($position_directory);\n";
			print "copy($log_file, catfile($position_directory, 'log.txt'));\n";
			print "&write_standard_config_file(catfile($position_directory,'analysis.cfg'),$exp_name)\n";
		} else {
			mkpath($position_directory);
			copy($log_file, catfile($position_directory, 'log.txt'));
			&write_standard_config_file(catfile($position_directory,'analysis.cfg'),$exp_name);
		}

		foreach my $field (keys %{$stage_position_sets{$position}}) {
			my $image_type = &find_image_type($field);
			
			my $no_spaces_field = $field;
			$no_spaces_field =~ s/ /_/g;

			my @temp = split("",@{$stage_position_sets{$position}{$field}});
			my $image_num_length = scalar(@temp);
			
			my $image_directory = catdir($position_directory, 'Images', $image_type);

			#we start the scan through each position field combination at one
			#because zero will always be empty because image numbering starts at
			#one
			foreach my $image_num (1 .. $#{$stage_position_sets{$position}{$field}}) {
				my $source_file = ${$stage_position_sets{$position}{$field}}[$image_num];
				
				my $suffix;
				if ($source_file =~ /.*\.(.*)/) {
					$suffix = $1;
				} else {
					die "Unable to find file suffix in: $source_file";
				}

				my $image_filename = $no_spaces_field . "_" . 
					sprintf('%0' . $image_num_length . 'd', $image_num) . 
					".$suffix";
				
				my $target_file = catfile($image_directory, $image_filename);
				$target_file =~ s/\.[tT][iI][fF][fF]?/\.png/;

				if ($image_num == 1) {
					if ($opt{debug}) {
						print "mkpath(dirname($target_file));\n";
					} else {
						mkpath(dirname($target_file));
					}
				}

				if ($opt{debug}) {
					print "convert $source_file $target_file;\n\n";
				} else {
					system("convert \"$source_file\" \"$target_file\" 2>/dev/null");
				}
			}
		}
	}
}

###############################################################################
# Functions
###############################################################################

sub find_image_type {
	my $field_name = $_[0];
	
	my $image_type = $field_translation{$field_name};
	if (defined $field_translation{$field_name}) {
		$image_type = $field_translation{$field_name};
	} else {
		die "Could not find appropriate way to translate image type for $field_name.";
	}
	
	return $image_type;
}

sub confirm_stage_hash {
	my %stage_hash = @_;
	
	#make sure the number of fields is the same in each position
	my $number_of_fields;
	foreach my $position (sort {$a<=>$b} keys %stage_hash) {
		if (defined $number_of_fields) {
			die "Different number of fields in stage position $position." 
				if $number_of_fields != scalar(keys %{$stage_hash{$position}});
		} else {
			$number_of_fields = scalar(keys %{$stage_hash{$position}});
		}
	}
	
	#make sure the number of images in each field in each position is the same
	foreach my $position (sort {$a<=>$b} keys %stage_hash) {
		my $image_num;
		foreach my $field (keys %{$stage_hash{$position}}) {
			if (defined $image_num) {
				die "Different number of images in each field in stage position $position." 
					if $image_num != scalar(@{$stage_hash{$position}{$field}});
			} else {
				$image_num = scalar(@{$stage_hash{$position}{$field}});
			}

			if (defined ${$stage_hash{$position}{$field}}[0]) {
				die "Expected the 0 index to be undefined, instead: ${$stage_hash{$position}{$field}}[0]";
			}
		}
	}
}

sub write_standard_config_file {
	my $file_name = shift;
	my $exp_name = shift;

	my $default_config_depth = &find_main_config_depth($file_name);
	my $config_prefix = "../" x $default_config_depth;
	

	open OUTPUT, ">$file_name";
	print OUTPUT 
"<<include " . $config_prefix . $relative_config_file . ">>

###############################################################################
# Required Experiment Parameters
###############################################################################

exp_name = $exp_name 

###############################################################################
# Experiment Specific Parameters
###############################################################################
";

	close OUTPUT;
}

sub determine_data_folder_name {
	my $config_file = $_[0];

	my %cfg = ParseConfig(
        -ConfigFile            => $config_file,
        -MergeDuplicateOptions => 1,
        -IncludeRelative       => 1,
	);
	
	return $cfg{data_folder};
}

sub find_main_config_depth {
	my $file_name = shift;

	my @split_file_parts = split(/\//,$file_name);
	
	while ($split_file_parts[0] eq "..") {
		@split_file_parts = @split_file_parts[1..$#split_file_parts];
	}
	
	return scalar(@split_file_parts) - 2;
}
