use strict;
use warnings;
use FindBin 1.51 qw( $RealBin );
use lib $RealBin;

use Getopt::Long;   # include the module for input

require data;


my $taxalist; #the input taxon list
my $opt = GetOptions( 'taxalist:s', \$taxalist); #set command line options
my $help;

if (!($opt && $taxalist) || $help) {#check for the required inputs
   print STDERR "\nExample usage:\n";
   print STDERR "\n$0 -taxalist=\"Danio_rerio Tetraodon_nigroviridis\" \n\n";
   print STDERR "Options:\n";
   print STDERR "    -taxalist = a list of subject species, use underscore to link genus and\n";
   print STDERR "                species name; the names should be separated by one space.\n\n";
   exit;
}


my @taxa = split(/\s/, $taxalist);

foreach my $taxon (@taxa) {

    my $infile1 = $taxon . "_R1.fq";
    my $infile2 = $taxon . "_R2.fq";
    my $newfile = $taxon . ".fas";
    my $newfile1 = $taxon . "_rmrep_R1.fq";
    my $newfile2 = $taxon . "_rmrep_R2.fq";
    my $index1 = $taxon . "_rmrep_R1.index";
    my $index2 = $taxon. "_rmrep_R2.index";

    my ($INFILE1, $INFILE2, $NEWFILE, $NEWFILE1, $NEWFILE2);

    open $INFILE1, "<$infile1" or die ("Cannot open $infile1 for reading ($!)");
    open $INFILE2, "<$infile2" or die ("Cannot open $infile2 for reading ($!)");
    open $NEWFILE, ">$newfile" or die ("Cannot open $newfile for writing ($!)");
    open $NEWFILE1, ">$newfile1" or die ("Cannot open $newfile for writing ($!)");
    open $NEWFILE2, ">$newfile2" or die ("Cannot open $newfile for writing ($!)");

    my ($id1, $id2);


    my %uniseq; #hash for storing unique sequence, hash key: 20bp first read + 20bp second read
    my ($key1, $key2, $key, $seq1, $seq2, $quality1, $quality2, $quality, $header1, $header2);

    while (my $line1 = readline ($INFILE1)) {

        chomp $line1;
        if ($line1 =~ /^@\S+\s+\S+/) { # if we find @
            if ($id1) { # if id exists we first need to write old information
                $key = "$key1" . "$key2";
                $quality = "$quality1" . "$quality2";
                if ($uniseq{$key}) {
                    my $qold = &data::averq($uniseq{$key}->{quality});
                    my $qnew = &data::averq($quality);
                    if ($qold < $qnew) {
                        $uniseq{$key} = { id => $id1,
                                          quality => $quality,
                                          seq1 => $seq1,
                                          seq2 => $seq2,
                                          header1 => $header1,
                                          header2 => $header2
                                        };
                    }
                }
                else {
                    if ($id1 eq $id2) {
                        $uniseq{$key} = { id => $id1,
                                          quality => $quality,
                                          seq1 => $seq1,
                                          seq2 => $seq2,
                                          header1 => $header1,
                                          header2 => $header2
                                        };
                    }
                     else {
                        die "$id1 and $id2 are not consistent\n";
                     }
                }
            }
            ($id1) = $line1 =~ /^@(\S+)/; #extract the new id
            $header1=$line1;
                $seq1 = readline ($INFILE1); #extract the sequence line as key for hash
                chomp $seq1;
                $key1 = substr($seq1, 0, 20); #get the first 20 nucleotides as key
                readline ($INFILE1);
                $quality1 = readline ($INFILE1);
                chomp $quality1;
        }

        #now we read INFILE2
        my $line2 = readline ($INFILE2);
        chomp $line2;
        if ($line2 =~ /^@\S+\s+\S+/) { # if we find @
            ($id2) = $line2 =~ /^@(\S+)/; #extract the new id
            $header2=$line2;
                $seq2 = readline ($INFILE2); #use the sequence line as key for hash
                chomp $seq2;
                $key2 = substr($seq2, 0, 20); #get the first 20 nucleotides as key
                readline ($INFILE2);
                $quality2 = readline ($INFILE2);
                chomp $quality2;
        }

    }


    #now make the fasta file without duplicats
    my @keys = keys (%uniseq);
    my @sortedkeys = sort (@keys);
    foreach my $key (@sortedkeys) {
            my $id = $uniseq{$key}->{id};
        my $header1 = $uniseq{$key}->{header1};
        my $header2 = $uniseq{$key}->{header2};
        my $seq1 = $uniseq{$key}->{seq1};
        my $seq2 = $uniseq{$key}->{seq2};

        my $q = $uniseq{$key}->{quality};
        my $q1 = substr($q, 0, length($seq1));
        my $q2 = substr($q, -length($seq2));

        print $NEWFILE ">$id\n$seq1\n$seq2\n";
        print $NEWFILE1 "\@$header1\n$seq1\n\+\n$q1\n";
        print $NEWFILE2 "\@$header2\n$seq2\n\+\n$q2\n";

    }


    close ($INFILE1) or die "Can't close the new file!!!";
    close ($INFILE2) or die "Can't close the new file!!!";
    close ($NEWFILE) or die "Can't close the new file!!!";
    close ($NEWFILE1) or die "Can't close the new file!!!";
    close ($NEWFILE2) or die "Can't close the new file!!!";

    &data::indexfasta($taxon); #make an index file for the new combined fasta file

    `makeblastdb -dbtype nucl -in $newfile -out $taxon`; #make a blast database for the new fasta file

    &data::indexfastq($newfile1, $index1); #index the read1 fastq

    &data::indexfastq($newfile2, $index2); #index the read2 fastq


}
