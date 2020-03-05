#!/usr/bin/perl -w
use Getopt::Long;
use File::Basename;
my $opt = GetOptions('query:s', \$query,
        'database:s', \$database);
if (!($opt && $query)) {#check for the required inputs
   print STDERR "\nExample usage:\n";
   print STDERR "\n$0 -query=\"Danio_rerio\" \n"; # query directory
   print STDERR "    -database = the name of blast database\n\n"; # database name
   exit;
}
`mkdir genebin`;
`mkdir blastout`;
my @fas=<$query/*.fas>;
foreach my $fas (@fas) {
        (my $gene)=basename($fas)=~/(.*)\.fas/;
        my @spec;
        my (%seq, %hash);
        my $spec;
        open "fil1", "$fas" or die "can not open $fas\n";
        while (<fil1>) {
                chomp;
                if (/>/) {
                        s/>//;
                        my @a=split;
                        $spec=$a[0];
                        push @spec, $spec;
                } else {
                        $seq{$spec}.=$_;
                }
        }
        `blastn -query $fas -db $database -out blastout/$gene.blast.out.txt -num_threads 32 -word_size 7 -gapopen 5 -gapextend 2 -penalty -1 -reward 1 -evalue 0.000001 -outfmt 6`;
        open "fil2", "blastout/$gene.blast.out.txt" or die "can not open $gene.blast.out.txt\n";
        while (<fil2>) {
                chomp;
                my @a=split;
                my $spec1=$a[0];
                my $gene1=$a[1];
                my $score=$a[-1];
                if ((defined $hash{$spec1} && $hash{$spec1}->{score}<=$a[-1]) || (!defined $hash{$spec1})) {
                        $hash{$spec1}={
                                gene => $gene1,
                                score => $score
                        };
                }
        }
        open fil3, ">genebin/$gene.fas" or die "can not open $gene.fas\n";
        foreach my $spec (@spec) {
                if ((defined $hash{$spec}) && ($hash{$spec}->{gene} eq $gene)) {
                        print fil3 ">$spec\n$seq{$spec}\n";
                }
        }
}
