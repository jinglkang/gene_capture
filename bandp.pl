#!/usr/bin/perl -w
#
###################################################################################################
#
# Name: bandp.pl
#
# This script is used to blast the query to non-redundant reads and parse the reads for different
# genes for each subject species.
#
# Datafile: fasta file query, local blast data base of the reads, indexed fastq files
#
# Input: the name of the query file
#
# Output: blastall output file and fastq reads file for different genes under each species
#
# Written by
#                 Chenhong Li
#                 Shanghai Ocean University, Shanghai, China.
#                 Created on Oct 2013
#                 Last modified on
#
###################################################################################################

use strict;
use warnings;


use Getopt::Long;   # include the module for input


my $query; #variable for store query species name
my $subject; #list of subject name
my $E = 0.0001; # E-value for blast, the default is 0.0001
my $help;

my $opt = GetOptions( 'query:s', \$query,
                      'subject:s', \$subject,
                      'E:f', \$E,
                      'help!', \$help); #set command line options


if (!($opt && $query) || $help) {#check for the required inputs
   print STDERR "\nExample usage:\n";
   print STDERR "\n$0 -query=\"Oreochromis_niloticus\" -subject=\"Rhyacichthys_aspro Odontobutis_potamophila\"\n\n";
   print STDERR "Options:\n";
   print STDERR "        -query = name of the query species \n";
   print STDERR "        -subject = a list of subject species, use underscore to link genus and\n";
   print STDERR "                     species name; taxa seperated by one space\n\n";
   print STDERR "For more optional parameters, see README.txt\n\n";
   exit;
}


my @subgenome = split (/\s/, $subject); #split the names of suber genomes
my @sortedsubgenome = sort @subgenome;


#call outside program blastall
my $queryfile = "$query.fas"; #get the name of blast query file
foreach my $sub (@sortedsubgenome) {
    my $databasefile = "$sub";
    my $blastout = "$query.$sub.blast.txt";
    print "Blasting, $query against $sub...\n";
    `blastn -query "$queryfile" -task blastn -db "$databasefile" -out "$blastout" -word_size 7 -gapopen 5 -gapextend 2 -penalty -1 -reward 1 -evalue $E -outfmt 6 -num_threads 32`;
    print "$query against $sub is done\n";



#now we need parse the blast results and extract original reads and bin them by query gene name for each species

#parse blast output files


    my $dir = $sub . "_results";
    `mkdir $dir`;

    #open the index file and fastq data file for retrieving late
    my ($SEQ_FILE1, $INDEX_FILE1);
    my $seqfile1 = $sub . "_rmrep_R1.fq";
    my $indexfile1 = $sub . "_rmrep_R1.index";
    open ($SEQ_FILE1, "<$seqfile1") or die "Cannot open $seqfile1 for reading ($!)";
    open ($INDEX_FILE1, "<$indexfile1") or die "Cannot open $indexfile1 for reading ($!)";
    #read the index1 file
    my %index1;
    my $id1;
    while (my $line = readline ($INDEX_FILE1)) {
        chomp $line;
        my ($id1, $position1) = split /\t/, $line;# the id in index file is longer than in embl file
        $index1{$id1} = $position1;
    }

    my ($SEQ_FILE2, $INDEX_FILE2);
    my $seqfile2 = $sub . "_rmrep_R2.fq";
    my $indexfile2 = $sub . "_rmrep_R2.index";
    open ($SEQ_FILE2, "<$seqfile2") or die "Cannot open $seqfile2 for reading ($!)";
    open ($INDEX_FILE2, "<$indexfile2") or die "Cannot open $indexfile2 for reading ($!)";
    #read the index2 file
    my %index2;
    my $id2;
    while (my $line = readline ($INDEX_FILE2)) {
        chomp $line;
        my ($id2, $position2) = split /\t/, $line;# the id in index file is longer than in embl file
        $index2{$id2} = $position2;
    }



    my $BLASTOUT;
    open ($BLASTOUT, "<$blastout") or die "Can't open the blastout file!!!"; #open blastout file for each subject
        my ($geneidlag, $hitidlag, $OUTFILE);

        my $line = readline ($BLASTOUT); #deal with the first line
        my ($geneid, $hitid) = $line =~ /^(\S+)\s+(\S+)/;
        my $outfile = $dir . "/" . "$geneid" . "." . "fq";
        open ($OUTFILE, ">$outfile") or die "Can't open the output file";
        $geneidlag = $geneid;
        $hitidlag = $hitid;

        while (my $line = readline ($BLASTOUT)) {

            ($geneid, $hitid) = $line =~ /^(\S+)\s+(\S+)/;

            if ($hitid ne $hitidlag) {

                #print $OUTFILE "$hitidlag\n";

                # and now we have superfast access to our sequences:
                seek $SEQ_FILE1, $index1{$hitidlag}, 0;  # let’s jump to our seq
                my $first = readline($SEQ_FILE1);  # first line is header, skip
                my $second = readline($SEQ_FILE1);      # second line is the whole sequence
                my $third = readline($SEQ_FILE1);       # third line is "+"
                my $fourth = readline($SEQ_FILE1);      # fourth line is the quality score
                my $seq = "$first" . "$second" . "$third" . "$fourth";
                print $OUTFILE "$seq";
                seek $SEQ_FILE2, $index2{$hitidlag}, 0;  # let’s jump to our seq
                $first = readline($SEQ_FILE2);  # first line is header, skip
                $second = readline($SEQ_FILE2); # second line is the whole sequence
                $third = readline($SEQ_FILE2);  # third line is "+"
                $fourth = readline($SEQ_FILE2); # fourth line is the quality score
                $seq = "$first" . "$second" . "$third" . "$fourth";
                print $OUTFILE "$seq";

            }

            if ($geneid ne $geneidlag){ #if $geneid is not equal to $geneidlag
                close ($OUTFILE) or die "Can't close the output file";
                $outfile = $dir . "/" . "$geneid" . "." . "fq";;
                open ($OUTFILE, ">$outfile") or die "Can't open the output file"; #open new output file
                $geneidlag = $geneid;
            }
            $hitidlag = $hitid;

        }

        #we need add the last line
        #print $OUTFILE "$hitidlag\n";


    close ($SEQ_FILE1) or die "Can't close the new file!!!";
    close ($INDEX_FILE1) or die "Can't close the index file!!!";
    close ($SEQ_FILE2) or die "Can't close the new file!!!";
    close ($INDEX_FILE2) or die "Can't close the index file!!!";

    #    `rm $sub*.* `   ;



}
