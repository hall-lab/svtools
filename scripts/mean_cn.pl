#!/usr/bin/env perl

use strict;
use POSIX qw(ceil floor);

my @samps;
my $min_freq = $ARGV[0];
my %dels     = (
    "0/0" => [],
    "0/1" => [],
    "1/1" => []
);
my %dups = (
    "0/0" => [],
    "0/1" => [],
    "1/1" => []
);

while ( my $line = <STDIN> ) {
    next if ( $line =~ /^##/ );
    chomp($line);
    if ( $line =~ /^#CHROM/ ) {
        @samps = split( /\s+/, $line );
    }
    else {
        my @line = split( /\s+/, $line );
        next if $line[7] =~ /SECONDARY/;
        my ($svtype) = $line[7] =~ /SVTYPE=([\w]+);/;
        if ( $svtype eq "DEL" || $svtype eq "DUP" ) {
            my @format = split( ":", $line[8] );
            my %index1;
            @index1{@format} = ( 0 .. $#format );
            my $gtind   = $index1{"GT"};
            my $cnind   = $index1{"CN"};
            my %gtcount = (
                "0/0" => 0,
                "0/1" => 0,
                "1/1" => 0
            );
            my %cnsum = (
                "0/0" => 0,
                "0/1" => 0,
                "1/1" => 0
            );
            for ( my $i = 9 ; $i < @samps ; $i++ ) {
                my @data = split( ":", $line[$i] );
                my $gt   = $data[$gtind];
                my $cn   = $data[$cnind];

                $cnsum{$gt} += $cn;
                $gtcount{$gt}++;
            }
            print "$line[2]\t$svtype\t";
            foreach my $key ( "0/0", "0/1", "1/1" ) {
                my $mean = $cnsum{$key} / ( $gtcount{$key} + 0.01 );
                print "$key\t$gtcount{$key}\t$mean\t";
                if ( $svtype eq "DEL" ) {
                    push @{ $dels{$key} }, $mean;
                }
                elsif ( $svtype eq "DUP" ) {
                    push @{ $dups{$key} }, $mean;
                }
            }
            print "\n";
        }
    }
}

foreach my $key ( "0/0", "0/1", "1/1" ) {
    my @temp = sort { $a <=> $b } @{ $dels{$key} };
    my $ll = scalar(@temp);
    print STDERR "DEL\t$key\t$temp[floor(0.1*$ll)]\t$temp[floor(0.9*$ll)]\t";
}
print STDERR "\n";

foreach my $key ( "0/0", "0/1", "1/1" ) {
    my @temp = sort { $a <=> $b } @{ $dups{$key} };
    my $ll = scalar(@temp);
    print STDERR "DUP\t$key\t$temp[floor(0.1*$ll)]\t$temp[floor(0.9*$ll)]\t";

}
print STDERR "\n";

