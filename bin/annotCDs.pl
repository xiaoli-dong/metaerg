
#!/usr/bin/env perl

use Getopt::Long;
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Benchmark;
use threads;
use threads::shared;
use List::Util qw(min max sum);
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::SearchIO;
require "util.pl";
require "search_parser.pl";

# Global variabies
my $starttime = localtime;
my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';
my $EXE = $FindBin::RealScript;
my $bin = "$FindBin::RealBin";

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# command line options
my(@Options, $outdir, $evalue, $cpus, $identity, $coverage, $hmm_cutoff, $hmm_evalue_cutoff, $DBDIR);
setOptions();
my $faa = shift @ARGV or err("Please supply a fasta format coding  sequence file on the command line.");

my @thrs = ();

push(@thrs, threads->new(\&cds_hmmerSearch_FOAM, $faa, "FOAM-hmm_rel1a.hmm", $outdir, $cpus, "--cut_tc"));
push(@thrs, threads->new(\&cds_diamondSearch,$faa, "uniprot_sprot", $outdir, $cpus, $evalue, 1));
push(@thrs, threads->new(\&cds_diamondSearch,$faa, "genomedb", $outdir, $cpus, $evalue, 0));
push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "Pfam-A.hmm", $outdir, $cpus,$hmm_cutoff, $hmm_evalue_cutoff));
push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "TIGRFAMs.hmm", $outdir, $cpus,$hmm_cutoff, $hmm_evalue_cutoff));

#push(@thrs, threads->new(\&cds_hmmerSearch_homemodel, $faa, "metabolic.hmm", $outdir, $cpus,$hmm_evalue_cutoff));
#push(@thrs, threads->new(\&cds_hmmerSearch_casgene, $faa, "casgenes.hmm", $outdir, $cpus,$hmm_evalue_cutoff));

push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "metabolic.hmm", $outdir, $cpus,"", $hmm_evalue_cutoff));
push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "casgenes.hmm", $outdir, $cpus,"", $hmm_evalue_cutoff));

foreach my $thr ( @thrs ) {
  $thr -> join();
}

my %seqHash = ();

#read all predicted features
open my $gff, "$outdir/features.gff" or die "could not open $outdir/features.gff to read, $!\n";
my @seqArray = ();
my $gffio = Bio::Tools::GFF->new(-fh =>$gff, -gff_version => 3);

while (my $f = $gffio->next_feature) {
    my $sid = $f->seq_id;
    push(@seqArray ,$sid) if not exists $seqHash{$sid};
    push (@{$seqHash{$sid}{FEATURE}}, $f);

}
close($gff);

parse_search_results($outdir, \%seqHash, $coverage, $identity, $hmm_evalue_cutoff);

#output features with db search annot
output_gff(\%seqHash, \@seqArray, $outdir);


sub output_gff{

    my ($seqHash, $seqArray, $outdir) = @_;
    my $gffver = 3;

    msg("Writing all feature gff file to $outdir");
    #open my $gff_fh, '>', "$outdir/$prefix.feature.gff";
    open my $gff_fh, '>', "$outdir/features.annot.gff";
    my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gffver);

    print $gff_fh "##gff-version $gffver\n";

    for my $sid (@seqArray) {
	for my $f ( sort { $a->start <=> $b->start } @{ $seqHash{$sid}{FEATURE} }) {
	    print $gff_fh $f->gff_string($gff_factory),"\n";
	}
    }
    close($gff_fh);
}



sub cds_diamondSearch{
    my ($faa, $dbname, $outdir, $cpus, $evalue, $mask) = @_;
    my $t0 = Benchmark->new;
    my $db = "$DBDIR/diamond/$dbname";
    my $bls_prefix = "$outdir/$dbname";
    my $blstable = "$bls_prefix\.blasttable";
    my $cmd = "diamond blastp --quiet -k 1 --masking $mask -p $cpus -q $faa -d $db -e $evalue --tmpdir /dev/shm -f 6 qseqid sseqid qlen qstart qend sstart send qframe pident bitscore evalue qcovhsp > $blstable 2> /dev/null";

    #msg("Will use diamond blastp search against $dbname:$cmd");
    if(! -e "$blstable"){
	runcmd("$cmd");
    }
    my $t1 = Benchmark->new;
    my $td = timediff($t1, $t0);
    msg("diamond blastp search using $dbname database took:" . timestr($td) . " to run\n");
}

sub cds_hmmerSearch_FOAM{
    my ($faa, $dbname, $outdir, $cpus, $hmm_cutoff) = @_;
    my $t0 = Benchmark->new;
    my $db = "$DBDIR/hmm/$dbname";
    my $hmmout = "$outdir/$dbname.tblout";
    my @thrs = ();
    push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "FOAM1.hmm", $outdir, $cpus, $hmm_cutoff));
    push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "FOAM2.hmm", $outdir, $cpus, $hmm_cutoff));
    push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "FOAM3.hmm", $outdir, $cpus, $hmm_cutoff));
    push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "FOAM4.hmm", $outdir, $cpus, $hmm_cutoff));
    push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "FOAM5.hmm", $outdir, $cpus, $hmm_cutoff));

    foreach my $thr ( @thrs ) {
	$thr -> join();
    }


    my @db_subset = ("FOAM1", "FOAM2","FOAM3","FOAM4","FOAM5");
    #my @db_subset = ("FOAM1", "FOAM2","FOAM3","FOAM4","FOAM5", "FOAM6","FOAM7","FOAM8", "FOAM9");
    my $cmd = "cat ";
    foreach my $db_subset_prefix (@db_subset){
	$cmd .= "$outdir/$db_subset_prefix.hmm.tblout ";
    }
    $cmd .= "> $hmmout";

    if(! -e $hmmout){
	runcmd($cmd);
    }
    my $t1 = Benchmark->new;
    my $td = timediff($t1, $t0);
    msg("hmmsearch using $dbname database took:" . timestr($td) . " to run\n");
}
sub cds_hmmerSearch{
    my ($faa, $dbname, $outdir, $cpus, $hmm_cutoff, $hmm_evalue_cutoff) = @_;
    my $t0 = Benchmark->new;
    my $db = "$DBDIR/hmm/$dbname";
    my $hmmprefix = "$outdir/$dbname";
    my $hmmout = "$hmmprefix.tblout";
    #my $cmd = "hmmsearch --notextw --acc $hmm_cutoff --cpu $cpus $db $faa  > \Q$hmmout\E 2> /dev/null";
    my $cmd = "";
    if($hmm_cutoff =~ /\S+/){
	$cmd = "hmmsearch --notextw --acc $hmm_cutoff --cpu $cpus --tblout \Q$hmmout\E $db $faa > /dev/null 2>&1";
    }
    else{
	$cmd = "hmmsearch --notextw --acc -E $hmm_evalue_cutoff --domE $hmm_evalue_cutoff --incE  $hmm_evalue_cutoff  --incdomE $hmm_evalue_cutoff --cpu $cpus --tblout \Q$hmmout\E $db $faa > /dev/null 2>\&1";
    }

    #msg("Will use hmmsearch against $dbname:$cmd");

    if(! -e $hmmout){
	runcmd($cmd);
    }
    my $t1 = Benchmark->new;
    my $td = timediff($t1, $t0);
    msg("hmmsearch using $dbname database took:" . timestr($td) . " to run\n");
}
#http://i.uestc.edu.cn/hmmcas/help.html


# Option setting routines

sub setOptions {
    use Getopt::Long;

    @Options = (
	'General:',
	{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
	{OPT=>"dbdir=s",  VAR=>\$DBDIR, DEFAULT=>"./db", DESC=>"metaerg searching database directory"},

	'Outputs:',
	{OPT=>"outdir=s",  VAR=>\$outdir, DEFAULT=>'', DESC=>"Output folder [auto]"},

	'Computation:',
	{OPT=>"cpus=i",  VAR=>\$cpus, DEFAULT=>8, DESC=>"Number of CPUs to use"},

	'diamond cutoff:',
	{OPT=>"evalue=f",  VAR=>\$evalue, DEFAULT=>1e-5, DESC=>"Similarity e-value cut-off"},
	{OPT=>"identity=f",  VAR=>\$identity, DEFAULT=>20, DESC=>"identity"},
	{OPT=>"coverage=f",  VAR=>\$coverage, DEFAULT=>70, DESC=>"coverage"},

	'hmmsearch cutoff:',
	{OPT=>"hmmcutoff=s",  VAR=>\$hmm_cutoff, DEFAULT=>"--cut_tc", DESC=>"hmm search trusted score threshold: [--cut_ga|--cut_nc|--cut_tc]"},
	{OPT=>"hmmevalue=f",  VAR=>\$hmm_evalue_cutoff, DEFAULT=>"1e-5", DESC=>"custom hmm db search threshold"}
	);

    (!@ARGV) && (usage());

    &GetOptions(map {$_->{OPT}, $_->{VAR}} grep { ref } @Options) || usage();

    # Now setup default values.
    foreach (@Options) {
	if (ref $_ && defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
	    ${$_->{VAR}} = $_->{DEFAULT};
	}
    }
}


#----------------------------------------------------------------------
sub usage {
    print STDERR
	"Name:\n  ", $EXE, " by $AUTHOR\n",
	"Synopsis:\n  coding sequence search against various database\n",
	"Usage:\n  $EXE [options] <fasta format protein coding sequence file>\n";
    foreach (@Options) {
	if (ref) {
	    my $def = defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	    $def = ($def ? ' (default OFF)' : '(default ON)') if $_->{OPT} =~ m/!$/;
	    my $opt = $_->{OPT};
	    $opt =~ s/!$//;
	    $opt =~ s/=s$/ [X]/;
	    $opt =~ s/=i$/ [N]/;
	    $opt =~ s/=f$/ [n.n]/;
	    printf STDERR "  --%-15s %s%s\n", $opt, $_->{DESC}, $def;
	}
	else {
	    print STDERR "$_\n";
	}
    }
    exit(1);
}


