# remove unconsistent sequences in paired fastq files
#!/usr/bin/perl -w
use Getopt::Long;
use File::Basename;
my $opt = GetOptions( 'taxalist:s', \$taxalist,
        'output:s', \$output); #set command line options
my @taxa = split(/\s/, $taxalist);
my $help;
if (!($opt && $taxalist) || $help) {#check for the required inputs
   print STDERR "\nExample usage:\n";
   print STDERR "\n$0 -taxalist=\"D13120286\" \n\n";
   print STDERR "Options:\n";
   print STDERR "       -taxalist = a list of fastaq files\n";
   print STDERR "       -output = output directory\n";
   exit;
}
foreach my $taxon (@taxa) {
        my (@header, @header1);
        my (%hash1, %hash2, %hash);
        my ($id, $seq, $quality);
    my $infile1 = $taxon . "_R1.fq"; # input paired fastq files
    my $infile2 = $taxon . "_R2.fq";
    my $infile1_1 = $taxon . "_R1.fq"; # output paired fastq files
    my $infile2_1 = $taxon . "_R2.fq";
    open "fil1", "<$infile1" or die "can't open $infile1";
    while (my $line1=readline(fil1)) {
        chomp $line1;
        if ($line1=~/^@\S+\s+\S+/) {
                ($id)=$line1=~/^@(\S+)/;
                push @header, $id; # push keys of %hash1 into @header
                $seq=readline (fil1);
                chomp $seq;
                readline (fil1);
                $quality=readline (fil1);
                chomp $quality;
                $hash1{$id}={
                        info => $line1,
                        quality => $quality,
                        seq => $seq
                };
        }
    }
    open "fil2", "<$infile2" or die "can't open $infile2";
    while (my $line1=readline(fil2)) {
        chomp $line1;
        if ($line1=~/^@\S+\s+\S+/) {
                ($id)=$line1=~/^@(\S+)/;
                push @header, $id; # push keys of %hash1 into @header
                $seq=readline (fil2);
                chomp $seq;
                readline (fil2);
                $quality=readline (fil2);
                chomp $quality;
                $hash2{$id}={
                        info => $line1,
                        quality => $quality,
                        seq => $seq
                };
        }
    }
    foreach my $header (@header) {
        $hash{$header}++;
        push @header1, $header if $hash{$header}==2; # if the key appeared in @header twice, then select it into @header1
    }
    unless (-f $output) {
        mkdir $output; # if there is no output directory, then mkdir it
    }
    foreach my $header (@header1) {
        open "fil1", ">>$output/$infile1_1" or die "cannot open $output/$infile1_1"; # make the output paired fastq files
        open "fil2", ">>$output/$infile2_1" or die "cannot open $output/$infile2_1";
        print fil1 "$hash1{$header}->{info}\n$hash1{$header}->{seq}\n+\n$hash1{$header}->{quality}\n" ; # print the info into new paired fastq files
        print fil2 "$hash2{$header}->{info}\n$hash2{$header}->{seq}\n+\n$hash2{$header}->{quality}\n" ; # print the info into new paired fastq files
    }
    %hash1=();
    %hash2=();
    %hash=();
    @header=();
}
