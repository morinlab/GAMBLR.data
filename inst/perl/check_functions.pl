#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Term::ANSIColor;
=begin
This script goes through all R code in the R directory and looks for functions that commonly have conflicts between tidyverse and another package and flags them
It's meant for developers to run to aid in avoiding bugs that can be hard to catch.
=cut 

chomp(my @scripts = `ls R/*.R`);
my @conflicts = ("select","filter","rename");

for my $script(@scripts){
    my @new_lines = flag_issues($script);
    print "\nAutomagically fix these issues in $script?\n\n";
    chomp(my $answer = <STDIN>);
    if($answer =~ /[yY]/){
        open OUT, ">$script";
        for(@new_lines){
            print OUT $_;
        }
        close OUT;
    }else{
        print "NOT CHANGING\n";
    }

}

sub flag_issues{
    my $script = shift;
    my $line_num = 0;
    my @lines;
    my $changed = 0;
    open S, $script or die "$!\n";
    while(my $line = <S>){
        $line_num++;
        for my $pattern (@conflicts){
            #ignore comments    
            my $comment = 0;
            if($line =~ /(.*)[^:]$pattern\(/){
                my $before = $1;
                if($before =~ /\#/){
                    $comment = 1;
                }
                if($comment){
                    print color('bold yellow');
                    print "CONSIDER DELETING COMMENTED LINE $line_num:\n$line";
                    print color('reset');
                    
                }
                else{
                    print color('bold red');
                    print "$script LINE $line_num matches $pattern\nOLD: $line";
                    print color('reset');
                    $line =~ s/$pattern\(/dplyr::$pattern\(/g;
                    print color('bold green');
                    print "NEW: $line";
                    
                    $changed++;
                    print color('reset');
                }

            }
            
        }
        push @lines, $line;
    }
    print "$changed lines of $line_num have suggested changes. Do you want to accept these? (Y/N)\n";
    
    
    return(@lines);
}