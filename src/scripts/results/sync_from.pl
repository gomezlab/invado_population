#!/usr/bin/perl -w

use Getopt::Long;

my %opt;
$opt{debug} = 0;
$opt{server} = "NOSERVER";
$opt{no_time} = 0;
$opt{repeat} = 1;
$opt{delay} = 0;
GetOptions(\%opt,"server=s", "debug|d", "no_progress|np", "no_time|nt", "repeat=s", "delay=s") or die;

my $progress_str = "--progress ";
if ($opt{no_progress}) {
    $progress_str = '';
}

my $time_str = "time ";
if ($opt{no_time}) {
    $time_str = '';
}

my $command = "$time_str" . "rsync $progress_str" . "-a $opt{server}:~/Documents/Projects/invadopodia/results/* ../../results/";
if ($opt{server} eq "NOSERVER" || $opt{debug}) {
    print "$command\n";
} else {
    while($opt{repeat} != 0) {
        system "$command";
        print "Done with sync $opt{repeat}\n\n";
        $opt{repeat}--;
        sleep $opt{delay}*60;
    }
}
