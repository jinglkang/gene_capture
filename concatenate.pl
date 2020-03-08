#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
# concatenate the sequences of all loci in each species
my $opt = GetOptions( 'clustal:s', \$clustal,
    'output:s', \$output);
if (!($opt && $clustal && $output) || $help) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 -clustal=\"clustal result directory\" -output=\"concatenate sequences file\"\n\n";
    exit;
}
my @fas=<$clustal/*.fas>;
my %seq;
my $spec;
foreach my $fas (@fas) {
	open "FAS", "$fas";
	while (<FAS>) {
		chomp;
		if (/>/) {
			s/>//;
			$spec=$_;
		} else {
			$seq{$spec}.=$_;
		}
	}
}
open "OUTPUT", ">$output";
foreach my $key (keys %seq) {
	print OUTPUT ">$key\n$seq{$key}\n";
}
