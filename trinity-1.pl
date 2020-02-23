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
        `find $results -name "*fq" | parallel -I% --max-args 1 Trinity --seqType fq --max_memory 30G --jaccard_clip --single % --full_cleanup --output %.Trinity`;
        my @fa=<*.fasta>;
        `rm *_map`;
        foreach my $fa (@fa) {
                (my $gene)=basename($fa)=~/(.*)\.fq/;
                (my $new_fa)=$gene.".fasta";
                `mv $fa $out/$new_fa`;
        }
}
