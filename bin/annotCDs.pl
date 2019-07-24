
#!/usr/bin/env perl

use Getopt::Long;
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin";
require "util.pl";
use Benchmark;
use threads;
use threads::shared;
use List::Util qw(min max sum);
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

# Global variabies
my $starttime = localtime;
my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';
my $VERSION = "0.1";
my $EXE = $FindBin::RealScript;
my $bin = "$FindBin::RealBin";
my $DBDIR = "$FindBin::RealBin/../db";
#my $programs = "$FindBin::RealBin/../programs";
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# command line options
my(@Options, $outdir, $evalue, $cpus, $identity, $coverage, $hmm_cutoff, $hmm_evalue_cutoff);
setOptions();
my $faa = shift @ARGV or err("Please supply a fasta format coding  sequence file on the command line.");

#cds_diamondSearch($faa, "genome-database", $outdir, $cpus, $evalue);
#cds_hmmerSearch($faa, "Pfam-A.hmm", $outdir, $cpus,$hmm_cutoff);

my @thrs = ();

push(@thrs, threads->new(\&cds_hmmerSearch_FOAM, $faa, "FOAM-hmm_rel1a.hmm", $outdir, $cpus, $hmm_cutoff));
push(@thrs, threads->new(\&cds_diamondSearch,$faa, "uniprot_sprot", $outdir, $cpus, $evalue, 1));
push(@thrs, threads->new(\&cds_diamondSearch,$faa, "genomedb", $outdir, $cpus, $evalue, 0));
push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "Pfam-A.hmm", $outdir, $cpus,$hmm_cutoff));
push(@thrs, threads->new(\&cds_hmmerSearch, $faa, "TIGRFAMs.hmm", $outdir, $cpus,$hmm_cutoff));
push(@thrs, threads->new(\&cds_hmmerSearch_homemodel, $faa, "metabolic.hmm", $outdir, $cpus,$hmm_evalue_cutoff));
push(@thrs, threads->new(\&cds_hmmerSearch_casgene, $faa, "casgenes.hmm", $outdir, $cpus,$hmm_evalue_cutoff));

foreach my $thr ( @thrs ) {
  $thr -> join();
}



sub cds_diamondSearch{
    my ($faa, $dbname, $outdir, $cpus, $evalue, $mask) = @_;
    my $t0 = Benchmark->new;
    my $db = "$DBDIR/diamond/$dbname";
    my $bls_prefix = "$outdir/$dbname";
    my $blstable = "$bls_prefix\.blasttable";
    my $cmd = "diamond blastp -k 1 --quiet -k 1 --masking $mask -p $cpus -q $faa -d $db -e $evalue --tmpdir /dev/shm -f 6 qseqid sseqid qlen qstart qend sstart send qframe pident bitscore evalue qcovhsp > $blstable 2> /dev/null";

    msg("Will use diamond blastp search against $dbname:$cmd");
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
    my $db = "$DBDIR/hmm/protein/$dbname";
    my $hmmout = "$outdir/$dbname.hmmer3";
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
	$cmd .= "$outdir/$db_subset_prefix.hmm.hmmer3 ";
    }
    $cmd .= ">> $hmmout";

    if(! -e $hmmout){
	runcmd($cmd);
    }
    my $t1 = Benchmark->new;
    my $td = timediff($t1, $t0);
    msg("hmmsearch using $dbname database took:" . timestr($td) . " to run\n");
}
sub cds_hmmerSearch{
    my ($faa, $dbname, $outdir, $cpus, $hmm_cutoff) = @_;
    my $t0 = Benchmark->new;
    my $db = "$DBDIR/hmm/protein/$dbname";
    my $hmmprefix = "$outdir/$dbname";
    my $hmmout = "$hmmprefix.hmmer3";
    my $cmd = "hmmsearch --notextw --acc $hmm_cutoff  --cpu $cpus $db $faa  > \Q$hmmout\E 2> /dev/null";
    msg("Will use hmmsearch against $dbname:$cmd");

    if(! -e $hmmout){
	runcmd($cmd);
    }
    my $t1 = Benchmark->new;
    my $td = timediff($t1, $t0);
    msg("hmmsearch using $dbname database took:" . timestr($td) . " to run\n");
}
sub cds_hmmerSearch_casgene{
    my ($faa, $dbname, $outdir, $cpus, $hmm_evalue_cutoff) = @_;
    my $t0 = Benchmark->new;
    my $db = "$DBDIR/hmm/protein/$dbname";
    my $hmmprefix = "$outdir/$dbname";
    my $hmmout = "$hmmprefix.hmmer3";
    my $cmd = "hmmsearch --notextw --acc -E $hmm_evalue_cutoff  --domT 25 --cpu $cpus $db $faa > \Q$hmmout\E 2> /dev/null";
    msg("Will use hmmsearch against $dbname database:$cmd");
    if(! -e $hmmout){
	runcmd($cmd);
    }
    my $t1 = Benchmark->new;
    my $td = timediff($t1, $t0);
    msg("hmmsearch using home made $dbname took:" . timestr($td) . " to run\n");
}
sub cds_hmmerSearch_homemodel{
    my ($faa, $dbname, $outdir, $cpus, $hmm_evalue_cutoff) = @_;
    my $t0 = Benchmark->new;
    my $db = "$DBDIR/hmm/protein/$dbname";
    my $hmmprefix = "$outdir/$dbname";
    my $hmmout = "$hmmprefix.hmmer3";
    my $cmd = "hmmsearch --notextw --acc -E $hmm_evalue_cutoff  --domT 25 --cpu $cpus $db $faa > \Q$hmmout\E 2> /dev/null";
    msg("Will use hmmsearch against $dbname database:$cmd");
    if(! -e $hmmout){
	runcmd($cmd);
    }
    my $t1 = Benchmark->new;
    my $td = timediff($t1, $t0);
    msg("hmmsearch using home made $dbname took:" . timestr($td) . " to run\n");
}

# Option setting routines

# Option setting routines

sub setOptions {
    use Getopt::Long;

    @Options = (
	'General:',
	{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
	{OPT=>"version", VAR=>\&version, DESC=>"Print version and exit"},

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
	"Name:\n  ", $EXE, " $VERSION by $AUTHOR\n",
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

#----------------------------------------------------------------------

sub version {
    print STDERR "$EXE $VERSION\n";
    exit;
}
