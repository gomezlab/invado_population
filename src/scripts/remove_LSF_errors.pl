#!/usr/bin/perl -w

use lib "../lib";
use lib "../lib/perl";
use Getopt::Long;
use File::Find::Rule;
use File::Spec::Functions;
use Data::Dumper;

use FA_job;
use Config::Adhesions qw(ParseConfig);

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "error_files=s") or die;

if (-e '/opt/lsf/bin/bjobs' && not $opt{lsf}) {
	die "LSF appears to be installed on this machine, don't you want to use it?" 
}	

if (not(exists($opt{cfg}) || exists($opt{error_files}))) {    
    die "Can't find cfg file specified on the command line";
}

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

$| = 1;
my $split_number = 45;

################################################################################
# Main
################################################################################

our @error_files;

if ($opt{error_files}) {
    @error_files = split(/,/,$opt{error_files});
    &remove_unimportant_errors(@error_files);
} else {
    &File::Find::find(\&ID_error_files, ($cfg{results_folder}));

    #split all the error file sets into the number of sets specified in
    #$split_number
    my @error_file_sets;
    my $set_pos = 1;
    while (@error_files) {
        my $temp = shift @error_files;
        push @{$error_file_sets[$set_pos - 1]}, $temp;
        $set_pos++;
        if ($set_pos > $split_number) {
            $set_pos = 1;
        }
    }
    
    #prepare the command for the LSF system
    my @commands;
    foreach my $set (@error_file_sets) {
        #$0 - the name of the program currently running, used to protect against
        #future file name changes
        push @commands, "$0 -lsf -error_files " . join(",",@{$set});
    }
    
    #these are needed for LSF to run correctly, but don't have anything to do
    #with the error files we have found using the above code, they will be used
    #to hold any error produced during the LSF run
    $opt{error_folder} = catdir('error_removal');

    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }

    &FA_job::send_general_lsf_program(\@commands,\%opt);

    exit(0);
}

################################################################################
# Functions
################################################################################

sub ID_error_files {
    if ($File::Find::name =~ /error.txt/) {
        push @error_files, $File::Find::name;
    }
}

sub remove_unimportant_errors {
    for (@_) {
        open INPUT, "$_" or die "$!: $_";
        my @errors = <INPUT>;
        close INPUT;

        my @cleaned_errors;
        foreach my $line (@errors) {
            if ($line =~ /Pending job threshold reached./) {
                next;
            }
            push @cleaned_errors, $line;
        }

        unlink $_;

        open OUTPUT, ">$_";
        print OUTPUT @cleaned_errors;
        close OUTPUT;
    }
    print "\n";
}
