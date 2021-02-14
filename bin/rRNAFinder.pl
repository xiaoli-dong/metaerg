#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Time::Piece;

use FindBin;
use Benchmark;
use List::Util qw[min max];
my $t0 = Benchmark->new;

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# global variables
my $EXE = $FindBin::RealScript;
my $VERSION = "1.1.0";
my $DESC = "genome/metagenome ribosomal RNA prediction and taxon assignments";
my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';
my $URL = 'https://sourceforge.net/projects/rrnafinder/';




# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# command line options
my(@Options, $quiet,$threads, $outdir, $evalue, $domain, $len, $identities, $coverage, $DBDIR);
setOptions();
my $bin = "$FindBin::RealBin";
$DBDIR ||= "$FindBin::RealBin/../db";

my $rna_hmmdir = "$DBDIR/hmm/rna";
if (! -e $outdir) {
    msg("Creating new output folder: $outdir");
    my $cmd = "mkdir -p \Q$outdir\E";
    system($cmd) >> 8 and  die "Could not execute cmd=$cmd, $!\n";
}

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# check all is well

msg("This is $EXE $VERSION");
msg("Written by $AUTHOR");
msg("Obtained from $URL");

my $fasta = shift @ARGV;
$fasta && -r $fasta or err("Usage: $EXE <file.fasta>");

my @result_lines = ();
my %read2hmmhit = ();

my @domains = ();
if($domain eq "meta"){
    push(@domains,("arc", "bac", "euk"));
}
elsif($domain =~ m/(arc|bac|euk)/){
    push(@domains, $domain);
}
else{
    err("--domain value can only be: arc|bac|euk|meta");
}


foreach my $kdom (@domains){

    my $hmmdb = "$rna_hmmdir/$kdom.hmm";
    
    $hmmdb && -r $hmmdb or err("Can't find database: $hmmdb");
    # run the external command
    my $search_output = "$outdir/$kdom.tblout";
    my $cmd = "nhmmer --cpu $threads -E $evalue -o /dev/null --tblout $search_output \Q$hmmdb\E \Q$fasta\E";
    
    msg("Scanning $fasta for $kdom rRNA genes using $hmmdb... please wait");
    msg("Command: $cmd");
    if(! -e $search_output){
	system($cmd) >> 8 and  die "Could not execute cmd=$cmd, $!\n";
    }
        
    # target name 0, accession 1, query name 2, accession 3, hmmfrom 4, hmm to 5, alifrom 6, ali to 7, envfrom 8, env to 9, sq len 10, strand 11, E-value 12, score 13, bias 14, description of target 15

    #for each location, only the best hit (lowest evalue and highest score) is reported
    open (OUT, $search_output) or die"Could not open $search_output to read, $!\n";
    foreach (<OUT>) {
	chomp;
	next if /^#/;    # comment line
	next if !/\S/; #empty line
	
	my @line = split(/\s+/, $_);
	err("bad line in nhmmer output - @line") unless defined $line[6] and $line[6] =~ m/^\d+$/;
	my($alifrom,$alito,$strand) = $line[6] < $line[7] ? ($line[6],$line[7],'+') : ($line[7],$line[6],'-');
	my($seqid, $qname, $evalue, $score) = ($line[0], $line[2], $line[12], $line[13]);
	
	next if($qname !~ /5S|5_8S/ && $score < 180);
			
	my $loc= $alifrom < $alito ? "$alifrom:$alito": "$alito:$alifrom";
	
	if(exists $read2hmmhit{$seqid}->{$loc}){
	    if($read2hmmhit{$seqid}->{$loc}->{evalue} > $evalue || $read2hmmhit{$seqid}->{$loc}->{score} < $score){
		$read2hmmhit{$seqid}->{$loc}->{kdom} = $kdom;
		$read2hmmhit{$seqid}->{$loc}->{evalue} = $evalue;
		$read2hmmhit{$seqid}->{$loc}->{score} = $score;
		$read2hmmhit{$seqid}->{$loc}->{alifrom} = $alifrom;
		$read2hmmhit{$seqid}->{$loc}->{alito} = $alito;
		$read2hmmhit{$seqid}->{$loc}->{strand} = $strand;
		$read2hmmhit{$seqid}->{$loc}->{qname} = $qname;
	    }
	    else{
		next;
	    }
	}
	#first time
	else{
	    $read2hmmhit{$seqid}->{$loc}->{kdom} = $kdom;
	    $read2hmmhit{$seqid}->{$loc}->{evalue} = $evalue;
	    $read2hmmhit{$seqid}->{$loc}->{score} = $score;
	    $read2hmmhit{$seqid}->{$loc}->{alifrom} = $alifrom;
	    $read2hmmhit{$seqid}->{$loc}->{alito} = $alito;
	    $read2hmmhit{$seqid}->{$loc}->{strand} = $strand;
	    $read2hmmhit{$seqid}->{$loc}->{qname} = $qname;
	}
    }
    close(OUT);
    #msg("Deleting unwanted file:", "$kdom.tblout");
    # unlink "$kdom.tblout";
}

#check overlap and delete the one with bigger evalue or smaller score when overlaps
foreach my $seqid (keys %read2hmmhit){
    my @dels = ();
    foreach my $loc1 (keys %{$read2hmmhit{$seqid}}){
	foreach my $loc2 (keys %{$read2hmmhit{$seqid}}){
	    next if $loc1 eq $loc2;
	    my ($s1, $e1) = split(/:/, $loc1);    
	    my($s2, $e2) = split(/:/, $loc2);
	    
	    if(is_overlapping($s1, $e1, $s2, $e2)){
		
		if($read2hmmhit{$seqid}->{$loc1}->{evalue} > $read2hmmhit{$seqid}->{$loc2}->{evalue} 
		   || $read2hmmhit{$seqid}->{$loc1}->{score} < $read2hmmhit{$seqid}->{$loc2}->{score}){
		    push(@dels, $loc1);
		}
		else{
		    push(@dels, $loc2);
		}
	    }
	}
    }
    foreach my $loc (@dels){
	delete $read2hmmhit{$seqid}->{$loc};
    }
    if(keys %{$read2hmmhit{$seqid}} == 0){
	delete $read2hmmhit{$seqid};
    }
}


my %rna2gff = ();
my %rnatcats = ();

foreach my $seqid (keys %read2hmmhit){
    my %rnaC = ();
    foreach my $loc (keys %{$read2hmmhit{$seqid}}){
	my $kdom = $read2hmmhit{$seqid}->{$loc}->{kdom};
	my $qname = $read2hmmhit{$seqid}->{$loc}->{qname};
	my $alifrom = $read2hmmhit{$seqid}->{$loc}->{alifrom};
	my $alito = $read2hmmhit{$seqid}->{$loc}->{alito};
	my $strand = $read2hmmhit{$seqid}->{$loc}->{strand};
	my $evalue = $read2hmmhit{$seqid}->{$loc}->{evalue};
	my $score = $read2hmmhit{$seqid}->{$loc}->{score};
	
	my $fragLen =  $alito-$alifrom;
	
	if($qname =~ m/(16S|18S|23S|28S)/ && $fragLen < $len){
	    msg("Skipping short rRNA: len=$fragLen $seqid\t$kdom\t$qname\t$alifrom\t$alito\t$evalue\t$strand\t$score\n");
	    next;
	}
	else{
	    $kdom =  $read2hmmhit{$seqid}->{$loc}->{kdom};
	    my ($cat) = $qname =~ /^\S+?\_(\S+)/;
	    $rnatcats{$cat}++;
	    $rnaC{$seqid}++;
	    my $ID = "$seqid\_" . $rnaC{$seqid};
	    my $gff_out_string = "$seqid\t$EXE\t$qname\t$alifrom\t$alito\t$evalue\t$strand\t\.\tID=$ID;Name=$qname;score=$score\n";
	    $rna2gff{$seqid}->{$ID} = $gff_out_string;
	}
    }
}


open (GFF, ">$outdir/rRNA.gff") or die "Could not open $outdir/rRNA.gff to write, $!\n";
my $gff_head = "##gff-version 3\n";
$gff_head .= "##$rna_hmmdir, $fasta\n";
print GFF $gff_head;

foreach my $sid (sort keys %rna2gff){
    foreach my $id (sort keys %{$rna2gff{$sid}}){
	print GFF $rna2gff{$sid}->{$id} if exists $rna2gff{$sid}->{$id};
    }
}
close(GFF);


gff2fasta("$outdir/rRNA.gff", $fasta, \%rnatcats);

if(! -e "$outdir/rRNA.tax.txt"){
    rna2taxon(\%rnatcats);
}

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print STDERR "rRNA prediction took:",timestr($td)," to run\n";

#Given two ranges [x1,x2], [y1,y2]
sub is_overlapping{
    my ($x1,$x2,$y1,$y2) = @_;
    return max($x1,$y1) <= min($x2,$y2)
}

sub rna2taxon{
    my ($rnatcats) = @_;
    
    my %taxfiles = ();
    foreach my $cat (sort keys %$rnatcats){
	my $input = "$cat.ffn";
	my $output = "$cat.tax.txt";
	$taxfiles{$cat} = "$outdir/$output";
	my $dbtype = "";
	if($cat =~ m/(16S|18S)/){
	    $dbtype = "ssu";
	}
	elsif($cat =~ m/(23S|28S|5S|5_8S)/ ){
	    $dbtype = "lsu";
	}
	else{
	    err("unknow rRNA type: $cat");
	    #next;
	}
	my $cmd = "$^X $bin/rna2taxon.pl --dbdir $DBDIR/blast --cpus $threads --dbtype $dbtype --evalue $evalue --identities $identities --coverage $coverage $outdir/$input > $outdir/$output";
	msg("Command: $cmd");
	system($cmd) >> 8 and  die "Could not execute cmd=$cmd, $!\n";
	
    }
    if(-e "$outdir/rRNA.tax.txt"){
	unlink("$outdir/rRNA.tax.txt");
    }
    foreach (keys %taxfiles){
	print STDERR "$_\n";
	my $catcmd = "cat $taxfiles{$_} >> $outdir/rRNA.tax.txt";
	
	msg("####################Command: $catcmd");
	system($catcmd) >> 8 and  die "Could not execute cmd=$catcmd, $!\n";
    }
    
}
sub gff2fasta{
    my ($gff, $fasta, $rnacats) = @_;
    
    my %rnas = ();
    open(GFF, $gff) or die "Could not open $gff to read, $!\n";
    open(FASTA, $fasta) or die "Could not open $fasta to read, $!\n";

    foreach (<GFF>){
	next if /^##/;
	my @a = split(/\s+/, $_);
	
       #9735487 rRNAFinder.pl   bac_16SrRNA        2       66      3.3e-12 -       .       ID=9735487_1;Name=bac_16SrRNA
	my $rnatype = $a[2];
	my ($id) =  $a[8] =~ /ID=(\S+?);/;
	$rnas{$a[0]}->{$rnatype}->{"$a[3]:$a[4]:$a[6]:$id"}=1;
    }
    close(GFF);

    my %file_handles = ();

    foreach my $cat (keys %$rnacats){
	msg("write to $outdir/$cat".".ffn");
	open ($file_handles{$cat}, ">$outdir/$cat".".ffn") or die "Could not open $outdir/$cat".".ffn to write, $!\n";
    }
    $/="\n>";
    while(<FASTA>){
	chomp;
	if(my ($sid,$other, $seq) =  /^>?(\S+)(.*?)\n(.*)/s){

	    if(exists $rnas{$sid}){
		$seq =~ s/\s+//g;
		foreach my $rnatype (keys %{$rnas{$sid}}){
		    my($cat) = $rnatype =~ /^\S+?\_(\S+)$/;
		    foreach my $loc (keys  %{$rnas{$sid}->{$rnatype}}){
			my ($start, $end, $strand, $id) = $loc =~ /(\d+)\:(\d+)\:(\S+)\:(\S+)/;
			my $length = $end-$start+1;
			if(($rnatype =~ /5S|5_8S/) || ($rnatype =~ m/(16S|18S|23S|28S)/ && $length >= $len)){
			    print {$file_handles{$cat}} ">$id type=$rnatype len=$length start=$start end=$end strand=$strand\n";
			    print  {$file_handles{$cat}} substr($seq, $start-1, $length), "\n";
			}

		    }
		}
	    }

	}

    }
}



#----------------------------------------------------------------------

sub msg {
  return if $quiet;
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print STDERR $line;
}

#----------------------------------------------------------------------

sub err {
  $quiet=0;
  msg(@_);
  exit(2);
}

#----------------------------------------------------------------------

sub version {
  print STDERR "$EXE $VERSION\n";
  exit;
}

# Option setting routines

sub setOptions {
  use Getopt::Long;
  @Options = (
      'Options:',
      {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
      {OPT=>"version", VAR=>\&version,           DESC=>"Print version and exit"},
      {OPT=>"quiet!",  VAR=>\$quiet, DEFAULT=>0, DESC=>"No screen output"},
      {OPT=>"dbdir=s",  VAR=>\$DBDIR, DEFAULT=>"", DESC=>"metaerg searching database directory"},
      
      {OPT=>"outdir=s",  VAR=>\$outdir, DEFAULT=>'.', DESC=>"Output folder [.]"},
      {OPT=>"threads=i",  VAR=>\$threads, DEFAULT=>8,  DESC=>"Number of threads/cores/CPUs to use"},
      
      {OPT=>"evalue=f",VAR=>\$evalue, DEFAULT=>1E-09, DESC=>"evalue cut-off for gene prediction and taxon assignment"},
      {OPT=>"domain=s",VAR=>\$domain, DEFAULT=>'meta', DESC=>"Search rRNA for rRNA in domains: [arc|bac|euk|meta]"},
      
      'rRNA gene prediction options',
      {OPT=>"length=s",  VAR=>\$len, DEFAULT=>180, DESC=>"length cut-off for 16/18/23/28S rRNA only"},
      
      'rRNA gene taxon assignment options:',
      {OPT=>"coverage=f",VAR=>\$coverage, DEFAULT=>80, DESC=>"The min percent of the query sequence ovelaping with the subject sequence"},
      {OPT=>"identities=s",VAR=>\$identities, DEFAULT=>'70,80,85,90,95,97,99', DESC=>"the percent identity of a match must exceed the given value of percent identity to be assigned at the given rank:domain70, phylum 80, class 85, order 90, family 95, genus 97, species 99"}
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

#
#----------------------------------------------------------------------
