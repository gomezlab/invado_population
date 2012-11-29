#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use lib "../lib";
use lib "../lib/perl";

use File::Basename;
use File::Spec::Functions;
use Cwd 'abs_path';
use Getopt::Long;
use File::Find::Rule;

use Config::Adhesions qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{extra} = "";
GetOptions(\%opt, "cfg|config=s", "debug|d", "program|p=s", "extra|e=s", 
                  "run_all_debug", "exp_filter=s", "no_parallel") or die;

die "Can't find cfg file specified on the command line" if not(exists $opt{cfg});
die "Can't find program to execute on the command line" if not(exists $opt{program});

print "Collecting Configuration\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################
chdir(dirname($opt{program}));

my $program_base = basename($opt{program});

my $debug_string = ($opt{debug}) ? "-d" : "";

my $cfg_suffix = basename($opt{cfg});
$cfg_suffix =~ s/.*\.(.*)/$1/;

my @config_files = File::Find::Rule->file()->name( "*.$cfg_suffix" )->in( ($cfg{data_folder}) );
@config_files = reverse sort @config_files;
if (exists($opt{exp_filter})) {
   @config_files = grep $_ =~ /$opt{exp_filter}/, @config_files;
}

my $parallel_return;
if ($opt{no_parallel}) {
	$parallel_return = 1;
} else {
	$parallel_return = system("which parallel");
}
my $ionice_return = system("ionice -c3 ls > /dev/null");

my @command_set;
foreach (@config_files) {
    next if /config\/default/;

    my $command = "nice -n 20 ./$program_base -cfg $_ $debug_string $opt{extra}";

	if ($ionice_return == 0) {
		$command = "ionice -c3 " . $command;
	}
	
    $command =~ s/"/\\"/g;
    push @command_set, "\"$command; echo $_;\"";
	if ($parallel_return != 0) {
		if ($opt{debug}) {
			print "$command\n";
		} else {
			print("$_\n");
			my $return_val = system("$command");
			if ($return_val != 0) {
				print("Caught Return Val: $return_val\n");
				last;
			}
		}
	}
}

if ($parallel_return == 0) {
	if ($opt{debug}) {
		my $parallel_cmd = "time parallel -u --nice 20 ::: \n\t" . join("\n\t", @command_set) . "\t\n";
		print $parallel_cmd;
	} else {
		while (@command_set) {
			my @sub_command_set;
			for (1..75) {
				if (@command_set) {
					push @sub_command_set, pop @command_set;
				}
			}
			my $parallel_cmd = "time parallel -u --nice 20 ::: " . join(" ", @sub_command_set);
			system($parallel_cmd);
		}
	}
}
