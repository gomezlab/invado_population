#!/usr/bin/perl

################################################################################
# Global Variables and Modules
################################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use threads;
use threads::shared;
use File::Spec::Functions;
use File::Basename;
use File::Find::Rule;
use Benchmark;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use POSIX;

use Config::ImageSet qw(ParseConfig);

my %opt;
$opt{debug} = 0;
$opt{lsf} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "exp_filter=s", "no_email",
	"sync=s") or die;

my $lsf_return = system("which bsub > /dev/null 2> /dev/null");

#Remember a return code of 0 means success
if ($lsf_return == 0 && not $opt{lsf}) {
	die "LSF appears to be installed on this machine, don't you want to use it?" 
}	

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

my $t1 = new Benchmark;
$|  = 1;

################################################################################
# Main
################################################################################

#This data structure acts as the overall control mechanism for which programs
#are run to perform the analysis. The first layer of the array are another set
#of arrays with sets of commands that can be executed simultaneously. The next
#layer holds all of those commands with the appropriate directory to execute the
#commands in.
my @overall_command_seq = (
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/find_median_images" ],
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/flat_field_correct_images" ],
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/find_exp_min_max" ],
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/find_cell_mask_experiment" ],
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/apply_bleaching_correction" ],
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/find_exp_min_max" ], 
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/find_cell_mask_properties_experiment" ],
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/track_cells_experiment" ],
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/find_cell_degrade_amount" ],
	[ ".", "./run_matlab_over_experiment.pl -script ../find_cell_features/measure_degradation_time_series" ],
	[ ".", "./run_matlab_over_experiment.pl -script ../analyze_cell_features/gather_tracking_results_experiment" ], 
	[ ".", "./run_matlab_over_experiment.pl -script ../analyze_cell_features/find_invading_cells_experiment" ], 
	[ ".", "./run_matlab_over_experiment.pl -script ../analyze_cell_features/build_overall_data_set.m" ], 
	[ ".", "./run_matlab_over_experiment.pl -script ../visualize_cell_features/build_single_cell_montage_experiment" ], 
	[ ".", "./run_matlab_over_experiment.pl -script ../visualize_cell_features/create_invader_visualization" ], 
	[ ".", "./run_matlab_over_experiment.pl -script ../visualize_cell_features/make_tracking_visualization" ], 
);

#some of the scripts only need to be run once for each experiment, this will
#rely on being able to find an experiment with "time_series_01" in its filename
my @run_only_once = qw(run_matlab_over_experiment);

my @skip_check = qw(find_median_images find_exp_min_max
	build_overall_data_set collect_montage_visualizations
	gather_tracking_results find_invading_cells track_cells
	build_single_cell_montage flat_field_correct_images
	find_cell_degrade_amount measure_degradation_time_series);

my $cfg_suffix = basename($opt{cfg});
$cfg_suffix =~ s/.*\.(.*)/$1/;

#config file processing
my @config_files = File::Find::Rule->file()->name( "*.$cfg_suffix" )->in( ($cfg{data_folder}) );
@config_files = sort @config_files;
if (exists($opt{exp_filter})) {
   @config_files = grep $_ =~ /$opt{exp_filter}/, @config_files;
}
die "No config files left after filtering." if (scalar(@config_files) == 0);

my @run_once_configs = grep $_ =~ /time_series_01/, @config_files;

#######################################
# Program Running
#######################################
my $starting_dir = getcwd;
for (@overall_command_seq) {
	my $command_start_bench = new Benchmark;
	my @command_seq = @{$_};
	@command_seq = map { $_->[0], $_->[1] . " -lsf" } @command_seq if $opt{lsf};
	print "Starting on $command_seq[1]\n";
	if (grep $command_seq[1] =~ /$_/, @run_only_once) {
		&execute_command_seq(\@command_seq, $starting_dir,\@run_once_configs);
	} else {
		&execute_command_seq(\@command_seq, $starting_dir);
	}

	#If debugging is not on, we want to wait till the current jobs finish
	#and then check the file complements of the experiments for completeness
	my $LSF_was_run = &wait_till_LSF_jobs_finish if ($opt{lsf} && not($opt{debug}));
	if (not($opt{debug}) && $LSF_was_run && not(grep $command_seq[1] =~ /$_/, @skip_check)) {
		print "Checking for all output files on command $command_seq[1]\n";
		my %exp_sets = &check_file_sets(\@config_files);

		for (1..2) {
			#if no experiments are left to retry, we break out and continue
			#to the next command set
			if (not(@{$exp_sets{retry}})) {
				last;
				next;
			}
			print "Retrying these experiments:\n".
				  join("\n", @{$exp_sets{retry}}) . "\n";
			&execute_command_seq(\@command_seq, $starting_dir, \@{$exp_sets{retry}});
			&wait_till_LSF_jobs_finish if ($opt{lsf});
			my %these_exp_sets = &check_file_sets(\@{$exp_sets{retry}});

			push @{$exp_sets{good}}, @{$these_exp_sets{good}};
			@{$exp_sets{retry}} = @{$these_exp_sets{retry}};
		}

		#check if there are any files left in the retry set, if so, clear
		#out the failed experiments from the next round
		if (@{$exp_sets{retry}}) {
			print "\nProblem with collecting full file complement on experiments after three tries:\n\t" .
				join("\n\t", @{$exp_sets{retry}}) . "\nRemoving them from the next run.\n\n";
			@config_files = @{$exp_sets{good}};  
		} else {
			print "Output file set complete, moving on.\n";
		}
	}

	if (not $opt{debug}) {
		my $command_end_bench = new Benchmark;
		my $td = timediff($command_end_bench, $command_start_bench);
		print "The command took:",timestr($td),"\n";
	}
	
	#we are done with the current command, so we print a few spacer lines to
	#indicate we have moved onto a new command
	print "\n\n";
}

if (not($opt{debug})) {
	my $t_bsub = new Benchmark;
	my $td = timediff($t_bsub, $t1);
	
	my $time_diff_str = "\"Took:" . timestr($td) . "\"";
	
	my $command;
	if ($opt{sync}) {
		$command = "results/sync_to.pl -server $opt{sync}";
		system($command);
	}
	
	$command = "bsub -J \"Job Finished: $opt{cfg}\" echo $time_diff_str";
	if (not $opt{no_email}) {
		system($command) if $opt{lsf};
	}
}

my $t2 = new Benchmark;
my $td = timediff($t2, $t1);
print "\nThe pipeline took:",timestr($td),"\n\n";

################################################################################
# Functions
################################################################################

sub execute_command_seq {
    my @command_seq  = @{ $_[0] };
    my $starting_dir = $_[1];
    my @these_config_files = @config_files;
    if (scalar(@_) > 2) {
        @these_config_files = @{$_[2]};
    }
	
	my $dir     = $command_seq[0];
	my $command = $command_seq[1];
	my $executed_scripts_count = 0;
	foreach my $cfg_file (@these_config_files) {
		my $config_command = "$command -cfg $cfg_file";
		chdir $dir;
		my $return_code = 0;
		$executed_scripts_count++;
		print "Done submitting: " if $executed_scripts_count == 1;

		if ($opt{debug}) {
			print $config_command, "\n";
		} else {
			$return_code = system($config_command);
			if ($return_code) {
				print "PROBLEM WITH: $config_command\n";
				print "RETURN CODE: $return_code\n";
			}
			if ($executed_scripts_count % ceil(scalar(@these_config_files)/10) == 0) {
				print sprintf('%.0f%% ',100*($executed_scripts_count/scalar(@these_config_files)));
			}
		}
		chdir $starting_dir;

		#if the return code was any number beside zero, indicating a problem
		#with the program exit, remove that config file from the run and
		#continue
		if ($return_code) {
			print "REMOVING: $cfg_file\n";
			@config_files = grep $cfg_file ne $_, @config_files;
		}
	}
	print "\n";
}

#######################################
# LSF
#######################################

sub wait_till_LSF_jobs_finish {
    #After each step of the pipeline, we want to wait till all the individual
    #jobs are completed, which will be checked three times
	my $total_checks = 1;
	my $first_check = &running_LSF_jobs;
	
	if (! $first_check) {
		print "No LSF jobs detected, jumping to next command.\n";
		return 0;
	} else {
		for (1 .. $total_checks) {
			print "LSF finished check number $_/$total_checks\n";
			my $sleep_time = 10;
			do {
				sleep($sleep_time);
			} while (&running_LSF_jobs);
		}
	}
	print "\n";
	return 1;
}

sub running_LSF_jobs {
    my $bjobs_command = "bjobs";
    if (defined $cfg{job_group}) {
        $bjobs_command .= " -g $cfg{job_group} 2>/dev/null";
    }

    my @lines = `$bjobs_command`;
    
    if (scalar(@lines) <= 1) {
        return 0;
    } else {
        return 1;
    }
}

sub check_file_sets {
    my @config_files = @{$_[0]};

    print "Out of " . scalar(@config_files) . ", # done:";
    my $number_configs_run = 0;
    
    my %exp_sets;
    #intialize these as the code that checks these matrices in the main program
    #assumes that a matrix will be present, even if it is empty
    @{$exp_sets{good}} = ();
    @{$exp_sets{retry}} = ();
    foreach my $this_config_file (@config_files) {
        my $return_code = system "./check_file_complement.pl -cfg $this_config_file";
        
        #if the return code is anything besides zero, add that config file back
        #to the exp_to_retry list
        if ($return_code) {
            push @{$exp_sets{retry}}, $this_config_file ;
        } else {
            push @{$exp_sets{good}}, $this_config_file ;
        }

        $number_configs_run++;
        if ($number_configs_run % ceil(scalar(@config_files)/10) == 0) {
            print " $number_configs_run";
        }
    }
    print "\n";

    return %exp_sets;
}
