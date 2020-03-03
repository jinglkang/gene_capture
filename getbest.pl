#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin 1.51 qw( $RealBin );
use lib $RealBin;
use smithwaterman; #module for pairwise alignment
use Getopt::Long;   #include the module for input

my $query; #variable for store query species name
my $subject; #list of subject name
my $flank = "n"; #print flank or n
my $help;

my $opt = GetOptions( 'query:s', \$query,
                      'subject:s', \$subject,
                      'f:s', \$flank,
                      'help!', \$help); #set command line options
                      

if (!($opt && $query) || $help) {#check for the required inputs
   print STDERR "\nExample usage:\n";
   print STDERR "$0 -query=\"Callorhinchus_milii\" -subject=\"Hydrolagus_affinis_1 Hydrolagus_affinis_2\"\n\n";
   print STDERR "Options:\n";
   print STDERR "        -query = name of the query species \n";
   print STDERR "        -subject = a list of subject species, use underscore to link genus and\n";
   print STDERR "                     species name; taxa seperated by one space\n";
   print STDERR "        -f = print flank or not \n";
   exit;
}

my @sub = split (/\s/, $subject); #split the names of subject taxon
chomp(@sub);
my @sortedsub = sort @sub;


#make dir for the output
my $resultdirnf = $query . ".resultnf";
my $resultnf = "result/$resultdirnf";
`mkdir $resultnf`;


my ($QUERY_FILE, $FASTA_FILENF, $FASTA_FILEF);
open $QUERY_FILE, "< $query.fas" or die ("Cannot open $query for reading ($!)");
while (my $line = readline ($QUERY_FILE)) {
    chomp $line;
    if ((my $gene) = $line =~ /^>(\S+)/) { # if we find >
        my $newnf = "result/$resultdirnf/$gene.fas"; #create files for each query genes
        open $FASTA_FILENF, ">$newnf" or die ("Cannot open $newnf for writing ($!)");
        my $seq1 = readline ($QUERY_FILE);
        chomp $seq1;
        
        #pairwise alignment: bait and assembled sequence from each subject 
        foreach my $subject (@sortedsub) {
            my $SUB_FILE;
            my $subgene = "trinity/$subject/" . $gene . ".fasta"; #find the gene file
            #print ">$subject\n";
            
            
            my $max_score;
            if (open $SUB_FILE, "<$subgene"){ #open the gene file of subject folder and do something, and if the gene file isn't exist, then escape 
                my (%score_seq);
                while (my $line = readline ($SUB_FILE)) {
                    chomp $line;
                    #print "$line\n";
                    if ((my $gene) = $line =~ /^>(\S+)/) {
                        my $seq2 = readline ($SUB_FILE);
                        chomp $seq2;
                        my $reseq2 = reverse $seq2;
                        $reseq2 =~ tr/ATCGatcg/TAGCtagc/;
                        my ($score, $newseq2) = &smithwaterman::smithwaterman($seq1, $seq2); #module for pairwise alignment
                        #print "$score, $newseq2\n";
                        $newseq2 =~ s/-//g; #replace the gap(-) with spacebar()
                        $score_seq{$score} = $newseq2;
                        ($score, $newseq2) = &smithwaterman::smithwaterman($seq1, $reseq2); #module for pairwise alignment
                        #print "$score, $newseq2\n";
                        $newseq2 =~ s/-//g; #replace the gap(-) with spacebar()
                        $score_seq{$score} = $newseq2;
                    }
             }   
                        my @score = sort {$a <=> $b} keys %score_seq;
                        #print "@score\n";
                        my $max_score = pop @score;
                        
                     print $FASTA_FILENF ">$subject\t$max_score\n$score_seq{$max_score}\n";
               
               %score_seq = ();
                              
               close ($SUB_FILE) or die "Can't close the fasta file!!!";
            }       
           }
         close ($FASTA_FILENF) or die "Can't close the fasta file!!!"; 
    }
   }
close ($QUERY_FILE) or die "Can't close the $query.fas file";
