#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use Time::Piece;
use Benchmark;
use Scalar::Util qw(openhandle);

my (@Options,$dbtype, $evalue, $cpus, $identities, $coverage,$length);
setOptions();
my $EXE = $FindBin::RealScript;
my $VERSION = "0.1";
my $DESC = "ribosomal RNA classification: associate ribosomal RNA to SILVA taxonomy";
my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';
my $OPSYS = $^O;

my $DBDIR = "$FindBin::RealBin/../db/blast";
my $ssu_db = "silva_SSURef_Nr99.fasta";
my $lsu_db = "silva_LSURef.fasta";
my $db = $dbtype eq "ssu" ? $ssu_db : $lsu_db;

my %level_cutoff = ();
my @arr = split(/,/, $identities);
if(@arr < 7){
    err("the number of the identity cutoff is less than 6, please redefine it" );
}
$level_cutoff{domain} = $arr[0];
$level_cutoff{phylum} = $arr[1];
$level_cutoff{class} = $arr[2];
$level_cutoff{order} = $arr[3];
$level_cutoff{family} = $arr[4];
$level_cutoff{genus} = $arr[5];
$level_cutoff{species} = $arr[6];

my $fasta = shift @ARGV;
$fasta && -r $fasta or err("Usage: $EXE <RNA fasta sequences file>");

my $t0 = Benchmark->new;

my $blastn_output = "$fasta.blastn.outfmt7.txt";
my $cmd = "blastn  -query $fasta -db $DBDIR/$db -dust no -num_threads $cpus -evalue $evalue -out $blastn_output  -outfmt \"7 qseqid qlen slen qstart qend sstart send length qcovhsp pident evalue bitscore stitle\" -max_target_seqs 5";

msg("$cmd\n");
if(! -e "$blastn_output"){
    runcmd($cmd);
}

open(BLAST_OUT,$blastn_output ) or die "Could not open $blastn_output to read, $!\n";

$/="\n# BLASTN 2.8.1+\n";
my $count = 0;
my %taxon_ass = ();

# Fields: query id 0, query length 1, subject length 2, q. start 3, q. end 4, s. start 5, s. end 6, alignment length 7, % query coverage per hsp 8, % identity 9, evalue 10, bit score 11, subject title 12
while(<BLAST_OUT>){
    # Query: NODE_278_length_54558_cov_93.7421_1 type=bac_16SrRNA len=180 start=1 end=758 strand=-
    my($qid, $rnatype) = $_ =~ /# Query:\s+?(\S+?)\s+type=(\S+)/;
    my @lines = split(/\n/, $_);
    my @hits = ();
    my $seq2taxon = "";
    foreach my $line (@lines){
	next if $line =~ /^#/;
	push(@hits, $line);
    }
    if (@hits == 0){
	$seq2taxon = "$qid\t$rnatype\tunknown\t[RNA_target=db:$db|nohit]";
	$taxon_ass{$qid} = $seq2taxon;
	next;
    }

    my $best_hit = $hits[0];
    my @best = split(/\t/, $best_hit);
    my $qlen = $best[1];
    my $cov = $best[8];
    my $iden = $best[9];
    
    if( ($rnatype !~ /5S|5_8S/) && ($qlen < $length || $cov < $coverage || $iden < $level_cutoff{phylum})){
	$seq2taxon = "$qid\t$rnatype\tunknown\t[RNA_target=db:$db|novalid]";
	$taxon_ass{$qid} = $seq2taxon;
    }
    else{
	$seq2taxon = getTaxon_by_bestHit($rnatype, $best_hit);
	$taxon_ass{$qid} = $seq2taxon;

    }
}

$/="\n";

foreach (keys %taxon_ass){

    print $taxon_ass{$_}, "\n";
}
my $t_end = Benchmark->new;
my $td = timediff($t_end, $t0);
msg("classifyRNA takes :",timestr($td)," to run\n");


sub getTaxon_by_bestHit{
    my ($rnatype,$best_hit) = @_;

# Fields: query id, query length, subject length, q. start, q. end, s. start, s. end, alignment length, % query coverage per hsp, % identity, evalue, bit score, subject title
    my @fields = split(/\t/, $best_hit);
    my $qid = $fields[0];
    my $qlen = $fields[1];
    my $iden = $fields[9];
    my $stitle = $fields[12];

    my ($ssid, $taxon) = $stitle =~ /(\S.*?)\s+\[(\S.*?)\]$/;
    my @taxon_terms = split(/;/, $taxon);

    my $target = "rRNA_target=db:$db|$ssid $fields[3] $fields[4] evalue:$fields[10] qcov:$fields[8] identity:$iden";
#the percent identity of a match must exceed the given value of percent identity to be assigned at the given rank: Species 99%, Genus 97%, Family 95%, Order 90%, Class 85%, Phylum 80%.
    #species
    my $taxon_assign = "$taxon";

    if($iden >= $level_cutoff{species} && $#taxon_terms > 6){
	$taxon_assign = join(";", @taxon_terms[0 .. 6]);
    }
    #genus
    elsif($iden >= $level_cutoff{genus} && $#taxon_terms > 5 ){
	$taxon_assign = join(";", @taxon_terms[0 .. 5]);
    }
    #family
    elsif($iden >= $level_cutoff{family} && $#taxon_terms > 4){
	$taxon_assign = join(";", @taxon_terms[0 .. 4]);
    }
    #order
    elsif($iden >= $level_cutoff{order} && $#taxon_terms > 3){
	$taxon_assign = join(";", @taxon_terms[0 .. 3]);
    }
    #class
    elsif($iden >= $level_cutoff{class} && $#taxon_terms > 2){
	$taxon_assign = join(";", @taxon_terms[0 .. 2]);
    }

    #phylum
    elsif($iden >= $level_cutoff{phylum} && $#taxon_terms > 1){
	$taxon_assign = join(";", @taxon_terms[0 .. 1]);
    }
    #domain
    elsif($iden >= $level_cutoff{domain} && $#taxon_terms >= 0){
	$taxon_assign = $taxon_terms[0];
    }
    return "$qid\t$rnatype\t$taxon_assign\t[$target]";
}

#the percent identity of a match must exceed the given value of percent identity to be assigned at the given rank: Species 99%, Genus 97%, Family 95%, Order 90%, Class 85%, Phylum 80%.
sub getTaxon_by_lca{
    my ($hits,$queryid, $qlen) = @_;

    my @taxons = ();
    my $index = 0;
    foreach (@$hits){
	# Fields: query id, query length, subject length, q. start, q. end, s. start, s. end, alignment length, % query coverage per hsp, % identity, evalue, bit score, subject title
	my @fields = split(/\t/, $_);

	my $cov = $fields[8];
	my $iden = $fields[9];
	next if $cov < $coverage;
	next if $iden < $level_cutoff{phylum};
	my ($taxon) = $fields[12] =~ /\[(\S.*?)\]$/;
	my @terms = split(/;/, $taxon);
	next if @terms == 0;
	@{$taxons[$index++]} = @terms;
    }
    my $com = lca(@taxons);
    return "$queryid\t$qlen\t" . join(";", @$com) . "\n";
}


sub lca {
    my @taxons = @_;

    my $size = @taxons;
    my $test = $taxons[0];

    #print STDERR "lca size=$size\n";
    #print STDERR "first*********=@$test\n";
    for (my $index = 0; $index <$size; $index++){
	#$test = $taxons[$index];
	#print STDERR "index=$index*********=@$test\n";
    }

    if (@taxons == 0) {return 0;} # Empty list? Not assigned=0

    @taxons=grep {defined} @taxons;  # screen out undefined values

    if (@taxons==0) {die("undefined values snuck into the taxid list; help?");} # how did undefined values get on the list?
    #print "  taxid3: ",(map {my @z=@$_;sprintf("%7s ",$z[$#z])} @taxons),"\n";

    if (@taxons==0) {return -1;} # No hits
    #print "  taxid4: ",(map {my @z=@$_;sprintf("%7s ",$z[$#z])} @taxons),"\n";

    # Find the lca
    my $n;
    for ($n=0;$n<100;$n++) {
	my @col;
	my $size = @taxons;

	for (my $m=0,@col=();$m<@taxons;$m++) { push @col,$taxons[$m][$n]; }
	if (!defined($col[0]) or !unanimous(@col)) {last;}
    }
    #if ($n<=0) {print STDERR "n=$n\n"; die("Impossible Error: Lca algorithm couldn't find any nodes in common!?  Data:\n".Dumper(\@_,\@taxons)." ");}

    my @slice = ();

    for(my $i = 0; $i < $n; $i++){
	push (@slice, $taxons[0][$i]);
    }

    #print STDERR "n=$n, common nodes=", $taxons[0][$n-1], "\n";
    #print STDERR "n*******=$n, common nodes=@slice\n";

    return \@slice;


    #return $taxons[0][$n-1]; #Return the taxid of the last common node
}
# Return 1 if the array is unanimous, 0 if not
sub unanimous {
  if (!defined $_[0]) {  # Unanimously 'undef'ined
    if (grep {defined} @_) {return 0;}
    return 1;
  }
  if (grep {!defined or $_ ne $_[0]} @_) {return 0;}
  return 1;
}



#----------------------------------------------------------------------

sub delfile {
    for my $file (@_) {

	msg("Deleting unwanted file:", $file);
	unlink $file or warn "Could not unlink $file: $!\n";;

    }
}
#----------------------------------------------------------------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print LOG $line if openhandle(\*LOG);
  print STDERR $line;
}

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

#species 99%, Genus 97%,Family 95%, Order 90%, Class 85%, Phylum 80%
sub setOptions {
  use Getopt::Long;
  @Options = (
      'Options:',
      {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
      {OPT=>"version", VAR=>\&version,           DESC=>"Print version and exit"},
      {OPT=>"cpus=i",  VAR=>\$cpus, DEFAULT=>8,  DESC=>"Number of threads/cores/CPUs to use"},
      {OPT=>"dbtype=s",VAR=>\$dbtype, DEFAULT=>'ssu', DESC=>"input sequence type:ssu|lsu, ssu incldues 16s/18s, lsu includes 23s/28s rRNA"},
      'Cutoffs:',
      {OPT=>"evalue=f",VAR=>\$evalue, DEFAULT=>1E-5, DESC=>"evalue cut-off"},
      {OPT=>"length=f",VAR=>\$length, DEFAULT=>200, DESC=>"rRNA gene length cutoff for taxon assignment"},
      {OPT=>"coverage=f",VAR=>\$coverage, DEFAULT=>80, DESC=>"The min query coverage: percent of the query sequence taht ovelaps the subject sequence"},
      {OPT=>"identities=s",VAR=>\$identities, DEFAULT=>'70,80,85,90,95,97,99', DESC=>"the percent identity of a match must exceed the given value of percent identity to be assigned at the given rank:domain 70, phylum 80, class 85, order 90, family 95, genus 97, species 99"},

      #{OPT=>"method=s",VAR=>\$method, DEFAULT=>'besthit', DESC=>"taxon assignment methods:besthit|lca|dominant"},
      #{OPT=>"toppercent=s",VAR=>\$toppercent, DEFAULT=>'25', DESC=>"the hits will be used for taxon assignment if the hits are within the distant of this value * best bitscore"},
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
  print STDERR "Synopsis:\n  $EXE $VERSION - $DESC\n";
  print STDERR "Author:\n  $AUTHOR\n";
  print STDERR "Usage:\n  $EXE [options] <contigs.fasta>\n";
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
sub decode {
  my $str = shift;
  $str =~ tr/+/ /;
  $str =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $str;
}
sub encode {
    my $str = shift;
    $str =~ s/([\t\n\r%&\=;,])/sprintf("%%%X",ord($1))/ge;
    return $str;
}
#----------------------------------------------------------------------

sub err {

  msg(@_);
  exit(2);
}
