# assemble acoording gene
#!/usr/bin/perl -w
use File::Basename;
use Getopt::Long;
my $opt = GetOptions( 'input:s', \$fq_dir,
        'output:s', \$output);
if (!($opt && $fq_dir && $output)) {
    print STDERR "\nExample usage:\n";
    print STDERR "\n$0 -input=\"fastq diretory\" /
    -output=\"trinity output direcotry\"\n\n";
    exit;
}
my @results=<$fq_dir/*_results>;
foreach my $result (@results) {
        (my $spec)=basename($result)=~/(.*)_results/;
        (my $out)="$output/$spec";
        mkdir $out;
        my @fq=<$result/*.fq>;
        foreach my $fq (@fq) {
                (my $gene)=basename($fq)=~/(.*)\.fq/;
                print "$fq\n";
                print "now is assembling $spec $gene ......\n";
                `Trinity --seqType fq --max_memory 30G --jaccard_clip --CPU 32 --single $fq --full_cleanup --output $gene.Trinity`;
                print "Trinity is over\n";
                my $file=$gene.".Trinity.Trinity.fasta";
                `mv $file $out/$gene.fasta`;
                `rm *_map`;
        }
}
