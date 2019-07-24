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
use Time::localtime;
use FileHandle;
#...........................................................................
# Global variabies
my $starttime = localtime;
my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';
my $VERSION = "0.1";
my $EXE = $FindBin::RealScript;
my $bin = "$FindBin::RealBin";
my $programs = "$FindBin::RealBin/../programs";
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# command line options
my(@Options, $prefix, $outdir, $force, $locustag, $gcode, $gtype, $minorflen, $sp, $tm, $cpus, $evalue, $mincontiglen);
setOptions();


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# prepare input, output filename, directory
$locustag ||= $EXE;
$prefix ||= $locustag.'_'.(localtime->mdy(''));
$outdir ||= $prefix;
if (!-d $outdir) {
    msg("Creating new output folder: $outdir\n");
    my $cmd = "mkdir -p \Q$outdir\E";
    (system $cmd) >> 8 and  err("Could not execute cmd=$cmd, $!\n");
}


my %seqHash = ();
my @seqArray = ();
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# prepare fasta, remove short contig, replace ambiguous with Ns
my $contigFile = shift @ARGV or err("Please supply a contig fasta file on the command line.");
my $fasta = "$outdir/$prefix.fna";
init($contigFile,$fasta, $mincontiglen, \%seqHash, \@seqArray);


# CRISPR
my $t0 = Benchmark->new;
my $mask = predict_CRISPRs($fasta, \%seqHash, $outdir);
maskSeqs($fasta, "$outdir/$prefix.crispr.masked.fna", $mask);
$fasta = "$outdir/$prefix.crispr.masked.fna";
my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
msg("predict_CRISPRs took:" . timestr($td) . " to run\n");

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# tRNA
$t0 = Benchmark->new;
$mask = predict_tRNA_aragorn($fasta,\%seqHash, $outdir);

maskSeqs($fasta, "$outdir/$prefix.tRNA.masked.fna", $mask);
$fasta = "$outdir/$prefix.tRNA.masked.fna";

$t1 = Benchmark->new;
$td = timediff($t1, $t0);
msg("predict_tRNA_aragorn took:" . timestr($td) . " to run\n");


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# rRNA: 16s, 18s, 23s, 28s, 5s, 5.8s
$t0 = Benchmark->new;
predict_rRNA($fasta, $gtype, \%seqHash, $cpus, $evalue, $outdir);
maskSeqs($fasta, "$outdir/$prefix.rRNA.masked.fna", $mask);
$fasta = "$outdir/$prefix.rRNA.masked.fna";

$t1 = Benchmark->new;
$td = timediff($t1, $t0);
msg("predict_rRNA took:" . timestr($td) . " to run\n");

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# CDS
$t0 = Benchmark->new;
my $num_cds = predict_CDs($fasta, \%seqHash,$prefix, $outdir, $gcode, $gtype);
$t1 = Benchmark->new;
$td = timediff($t1, $t0);
msg("predict_CDs took:" . timestr($td) . " to run\n");

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# Connect features to their parent sequences
msg("Connecting features back to sequences");
for my $sid (@seqArray) {
    for my $f (@{$seqHash{$sid}{FEATURE}}) {
	$f->attach_seq( $seqHash{$sid}{DNA} );
    }
    
}

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# predicts the presence and location of signal peptide cleavage sites in amino acid sequences from different organisms:
#Gram-positive prokaryotes, Gram-negative prokaryotes, and eukaryotes.
#The method incorporates a prediction of cleavage sites and a signal peptide/non-signal peptide prediction
#based on a combination of several artificial neural networks.

my $cdsoutfn = "$outdir/signalp.faa";

if($sp){
    $t0 = Benchmark->new;
    predict_signal_peptide("$outdir/cds.faa", \%seqHash, \@seqArray);
    $t1 = Benchmark->new;
    $td = timediff($t1, $t0);
    msg("predict_signal_peptide took:" . timestr($td) . " to run\n");
}

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# predict transmembrane helices based on a hidden Markov model on predicted CDs
$cdsoutfn = "$outdir/tmhmm.faa";
if($tm){
    $t0 = Benchmark->new;
    predict_TM("$outdir/cds.faa", \%seqHash, \@seqArray);
    $t1 = Benchmark->new;
    $td = timediff($t1, $t0);
    msg("predict_TM took:" . timestr($td) . " to run\n");
    
}

output_gff(\%seqHash, \@seqArray, $prefix, $outdir);
#output_fasta( \%seqHash, \@seqArray, $prefix, $outdir);

#----------------------------------------------------------------------
sub predict_tRNA_aragorn{
    
    my ($fasta, $seq, $outdir) = @_;
    my %mask = ();
    msg("Predicting tRNAs");
    
    #default is metagenome mode
    my $cmd = "aragorn -l -t -gc11  $fasta -w -o $outdir/tRNA.temp";
    
    if(! -e "$outdir/tRNA.temp"){		  
	runcmd($cmd);
    }				  
    open TRNA, "$outdir/tRNA.temp";
    my $num_trna=0;
    
#open TRNA, '-|', $cmd;
    
    #>C1997739
    #0 genes found
    #>C1997905
    #2 genes found
    #1   tRNA-Tyr               [170962,171044]      35      (gta)
    #2   tRNA-???               [172294,172367]      35      (aa)
    #1   tRNA-Leu                       [-3,82]      35      (caa)
    #1   tmRNA*                       c[-31,312]      208,252 ANDNAQTGAVALAA*

    #>end    957850 sequences 11875 tRNA genes
    $/ = "\n>";
    while (<TRNA>) {
	chomp;
	last if /^>?ebd/;
	if(my ($fastaName,$found) =  /^>?(\S+).*?\n(.*)/s){
	    next if $found =~ /0 genes found/;
	    
	    my @hits = split(/\n/, $found);
	    shift @hits;
	    my $sid = $fastaName;
	    my $strand = +1;
	    foreach my $tRNA (@hits){
		my @line = split(/\s+/, $tRNA);
		my $tcount = $line[0];
		$strand = -1 if $tRNA =~ /c\[/;
		
		my ($start, $end) = $line[2] =~ /\[(-?\d+?),(-?\d+)\]/;
		# correct strange coordinates in -l mode
		$start = max( $start, 1 );
		$end = min( $end, $seq->{$sid}{DNA}->length );
		
		my $tRNAType = $line[1];
		$tRNAType =~ s/-/_/g;
		my ($antiCodon) = $line[4] =~ /\((\S+?)\)/;
		
		$num_trna++;
		my $product = $tRNAType."_" . $antiCodon;
		my $tool = "aragorn";
		push @{$seq->{$sid}{FEATURE}},Bio::SeqFeature::Generic->new(
		    -primary    => "tRNA",
		    -seq_id     => $sid,
		    -source     => $tool,
		    -start      => $start,
		    -end        => $end,
		    -strand     => $strand,
		    #-score      => $cove_score,
		    -frame      => ".",
		    -tag        => {
			'Name' => $product,
			'ID' => "$sid\_tRNA\_$tcount"
		    }
		    );
		$mask{$sid}->{"$start:$end"} = 1;
	    }
	}
    }
    msg("Found $num_trna tRNAs");
    $/ = "\n";
    close(TRNA);
    return \%mask;
    #return $seq;
    
}

#----------------------------------------------------------------------
sub predict_rRNA{

    my ($fasta, $gtype, $seq, $cpus, $evalue, $outdir) = @_;
    my %mask = ();
    
    #9764223 nhmmer_rRNA.pl  ssu_rRNA        1098    1150    1.5e-10 +       .       Name=euk_18SrRNA
    msg("Predicting Ribosomal RNAs");
    my $cmd = "$bin/rRNAFinder.pl --threads $cpus --evalue $evalue --domain meta --outdir $outdir $fasta";
    
    #my $cmd = "rRNAFinder.pl --threads $cpus --evalue $evalue --domain $gtype --outdir $outdir $fasta";
    
    if(! -e "$outdir/rRNA.gff"){
	runcmd($cmd);
    }
    #NODE_161294_length_1378_cov_32.6674_1	16SrRNA	Bacteria;Planctomycetes;vadinHA49	[rRNA_target=db:silva_SSURef_Nr99.fasta|EF632951.1.1540 1 685 evalue:0.0 qcov:98 identity:92.453]
    my $num_rrna=0;
    
    if(! -e "$outdir/rRNA.tax.txt"){
	msg("Found $num_rrna rRNAs");
	return \%mask;

    }
    
    open(TAX, "$outdir/rRNA.tax.txt") or die "Could not open $outdir/rRNA.tax.txt to read, !$\n";
    my %rna2taxon = ();
    while(<TAX>){
	chomp;
	next if!/\S/;
	my @l = split(/\t/, $_);
	my $type = $l[1];
	my $taxon = $l[2];
	my ($target) = $l[3] =~ /\[(\S.*?)\]/;
	$rna2taxon{$type}->{$l[0]}->{taxon} = $taxon;
	$rna2taxon{$type}->{$l[0]}->{target} = $target;

	   
    }
    
    open my $NHMMER_RRNA, "$outdir/rRNA.gff";
        
    my $gff = Bio::Tools::GFF->new(-fh => $NHMMER_RRNA, -gff_version => 3);
    while (my $f = $gff->next_feature) {
	
	my $primary_tag = $f->primary_tag;
	my $feature_id = ($f->get_tag_values("ID"))[0];
	#my ($tmp) = $primary_tag =~ /^\w+?_(\w+)/;
	my $tmp = $primary_tag;
	if(exists $rna2taxon{$tmp}->{$feature_id}){
	    $f->add_tag_value("rRNA_taxon", $rna2taxon{$tmp}->{$feature_id}->{taxon});
	    $f->add_tag_value('rRNA_target', $rna2taxon{$tmp}->{$feature_id}->{target});
	}
	$f->score(".");
	my $sid = $f->seq_id;
	my $start = $f->start;
	my $end = $f->end;
	
	push @{$seq->{$sid}{FEATURE}}, $f;
	$mask{$sid}->{"$start:$end"} = 1;
	
	$num_rrna++;
    }
    msg("Found $num_rrna rRNAs");
    close($NHMMER_RRNA);
    return \%mask;
}


#----------------------------------------------------------------------
#MinCED is a program to find Clustered Regularly Interspaced Short Palindromic
#Repeats (CRISPRs) in full genomes or environmental datasets such as metagenomes,
#in which sequence size can be anywhere from 100 to 800 bp.
sub predict_CRISPRs{

    my($fasta, $seq, $outdir) = @_;
    my %mask = ();
    
    msg("Searching for CRISPR repeats");
    my $num_crispr=0;
    my $cmd = "minced -gffFull $fasta $outdir/crisprs.temp";
        
    if(! -e "$outdir/crisprs.temp"){
    	 runcmd($cmd);
    }
    open my $MINCED, "$outdir/crisprs.temp";
    #open my $MINCED, '-|', $cmd;
    my $gff = Bio::Tools::GFF->new(-fh => $MINCED, -gff_version => 3);
    #Sequence tool CRISPR crispr.start crispr.end numRepeat_CRISPR 
    #NODE_883_length_34292_cov_11.1505	minced:0.3.2	repeat_region	33687	33944	4	.	.	ID=CRISPR3;bin=91;rpt_family=CRISPR;rpt_type=direct;rpt_unit_seq=GTCGCACTGGGCTTCTAAAGCCCATGAGGATTGAAAC
    #NODE_883_length_34292_cov_11.1505	minced:0.3.2	repeat_unit	33687	33723	1	.	.	ID=DR.CRISPR3.1;Parent=CRISPR3;bin=91
    while (my $f = $gff->next_feature) {
	next if $f->primary_tag ne "repeat_region";
	my $sid = $f->seq_id;
	my $start = $f->start;
	my $end = $f->end;
	my $spacer_count = $f->score;
	$f->add_tag_value("note", "CRISPR with $spacer_count repeat units");
	$f->score(".");
	my $id = TAG($f, "ID");
	
	$mask{$sid}->{"$start:$end"} = 1;
		
	push @{$seq->{$sid}{FEATURE}}, $f;
	$num_crispr++;
	
    }
    msg("Found $num_crispr CRISPRs");
    close($MINCED);
    return \%mask;

}


#----------------------------------------------------------------------
sub predict_CDs{

    my($fasta, $seqHash,$prefix, $outdir, $gcode, $gtype) = @_;
    msg("Predicting coding sequences");
    #my $cmd = "prodigal -p meta -m -f gff -q -i $fasta ";
    my $proc = "meta";
    if($gtype ne "meta"){
	$proc = "single";
    }
    
    my $cmd = "prodigal -g $gcode -p $proc -m -f gff -q -i $fasta -a $outdir/cds.faa.temp.1 -d $outdir/cds.ffn.temp.1 -o $outdir/cds.gff.temp.1";
    
    if(! -e "$outdir/cds.gff.temp.1" || ! -e "$outdir/cds.faa.temp.1" || ! -e "$outdir/cds.ffn.temp.1"){		  
	runcmd($cmd);
    }
    open my $PRODIGAL, "$outdir/cds.gff.temp.1"; 
    #open my $PRODIGAL, '-|', $cmd;
    my $gff = Bio::Tools::GFF->new(-fh => $PRODIGAL, -gff_version => 3);
    $num_cds = 0;
    
    #scaffold7       Prodigal_v2.6.1 CDS     39      407     29.7    -       0       ID=2_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.415;conf=99.89;score=29.70;cscore=18.64;sscore=11.06;rscore=6.02;uscore=1.46;tscore=3.58;
    #scaffold7       Prodigal_v2.6.1 CDS     503     904     27.4    -       0       ID=2_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.463;conf=99.79;score=26.78;cscore=19.00;sscore=7.78;rscore=1.61;uscore=3.24;tscore=3.58;
    my $num_lt=0;
    my %excludes = ();
    while (my $f = $gff->next_feature) {
	foreach my $tag (my @tags = $f->get_all_tags()){
	    $f->remove_tag($tag) if ($tag ne "ID" &&  $tag ne "score");
	}
	my $sid = $f->seq_id;
	my $start = $f->start;
	my $end = $f->end;
	my $strand = $f->strand;
	#$f->add_tag_value('Parent', $sid);
	
	if(abs($end - $start) < $minorflen){
	    $excludes{TAG($f, "ID")} = 1;
	    msg("Excluding CDS too short $sid ". TAG($f, "ID"));
	}
	else{
	    my ($newid) = TAG($f, "ID") =~ /^\d+?_(\d+)/;
	    foreach my $tag (my @tags = $f->get_all_tags()){
		$f->remove_tag($tag); #if ($tag ne "ID" &&  $tag ne "score");
	    }
	    $f->add_tag_value("ID", "$sid\_cds\_$newid");
	    #print STDERR "newid = $sid\_$newid\n";
	    push @{$seqHash->{$f->seq_id}{FEATURE}}, $f;
	    $num_cds++;
	}
    }
    msg("Found $num_cds CDS");

    #>NODE_17_6 # 5936 # 6664 # 1 # ID=1_6;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.571
    #exclue the genes which are too short
    open AA_TEMP, "$outdir/cds.faa.temp.1"  or err("Can't open $outdir/cds.faa.temp.1");
    open NT_TEMP,  "$outdir/cds.ffn.temp.1" or err("Can't open $outdir/cds.ffn.1.temp");
    
    #the outfiles have already excluded some sequences
    open AA, ">$outdir/cds.faa"  or err("Can't open $outdir/cds.faa"); 
    open NT,  ">$outdir/cds.ffn" or err("Can't open $outdir/cds.ffn");
    
    $/ = "\n>";
    
    while(<AA_TEMP>){
	chomp;
	if(my ($id, $header,$seq) =  /^>?(\S+?)\s+(\S.*?)\n(.*)/s){
	    my @header_items = split(/\s+/, $header);
	    my ($ID) = $header_items[7] =~ /^ID=(\S+?);/;
	   
	    $id =~ s/^(\S+)\_(\d+)/$1\_cds\_$2/;

	    $seq =~ s/\*$//;
	    print AA ">$id $header\n$seq\n" if(not exists $excludes{$ID});
	}
    }
    while(<NT_TEMP>){
	chomp;
	if(my ($id, $header,$seq) =  /^>?(\S+?)\s+(\S.*?)\n(.*)/s){
	    
	    my @header_items = split(/\s+/, $header);
	    my ($ID) = $header_items[7] =~ /^ID=(\S+?);/;
	    $id =~ s/^(\S+)\_(\d+)/$1\_cds\_$2/;
	    
	    print NT ">$id $header\n$seq\n" if(not exists $excludes{$ID});
	}
    }
    close(AA_TEMP);
    close(NT_TEMP);
    
    close(AA);
    close(NT);
    
    $/ = "\n";
    #delfile("$outdir/cds.faa.temp", "$outdir/cds.ffn.temp");
    return $num_cds;
}


#----------------------------------------------------------------------
sub predict_signal_peptide{
    my ($fasta, $seqHash, $seqArray) = @_;
    msg("Looking for signal peptides at start of predicted proteins in meta mode");
    
    #open my $cdsoutfh, '>', $fasta;
    #my $cdsout = Bio::SeqIO->new(-fh=>$cdsoutfh, -format=>'fasta');
    my %cds;
    my $signalp = "signalp";
    for my $sid (@$seqArray) {
	for my $f (@{ $seqHash->{$sid}{FEATURE} }) {
	    next unless $f->primary_tag eq "CDS"; 
	    my ($id) = TAG($f, "ID");
	    #print STDERR "------id=$id\n";
	    $cds{$id} = $f;
	}
    }
    
    if(! -e "$fasta.signalp.txt.temp"){
	my $cmd1 = "$signalp -t gram- -f summary $fasta > $fasta.signalp.gramn.txt.temp;";
	my $cmd2 = "$signalp -t gram+ -f summary $fasta > $fasta.signalp.gramp.txt.temp;";
	my $cmd3 = "$signalp -t euk -f summary $fasta > $fasta.signalp.euk.txt.temp;";

	my $thr1 = threads->new(\&runcmd, $cmd1);
	my $thr2 = threads->new(\&runcmd, $cmd2);
	my $thr3 = threads->new(\&runcmd, $cmd3);
	$thr1->join();
	$thr2->join();
	$thr3->join();
	
	msg("Running: cat $fasta.signalp.euk.txt.temp $fasta.signalp.gramp.txt.temp $fasta.signalp.gramn.txt.temp > $fasta.signalp.txt.temp");
	system("cat $fasta.signalp.euk.txt.temp $fasta.signalp.gramp.txt.temp $fasta.signalp.gramn.txt.temp > $fasta.signalp.txt.temp");
    }
    my %peptides = ();
    my $tool = "SignalP";
    open(SIGNALP, "$fasta.signalp.txt.temp") or die "Could not open $fasta.signalp.txt to read, $!\n";
    #Name=METAANNOT_00002    SP='YES' Cleavage site between pos. 22 and 23: AHA-QT D=0.799 D-cutoff=0.570 Networks=SignalP-noTM
   
    my $num_sigpep = 0;
    $/ ="\n# Measure";
    my $sp_count = 0;
    while (<SIGNALP>) {
	next if !/SP=\'YES\'/;
	next if !/Name/;
	if( my ($spstart, $spend,$id, $dscore, $dcutoff, $networks) = $_ =~ /\n\s+D\s+(\d+)-(\d+)\s+.*?\nName=(\S+)\s+.*?D=(\S+)\s+D-cutoff=(\S+)\s+Networks=(\S+)/s){
	    my $parent = $cds{$id};
	    if(not defined $parent){
		print STDERR "************************************* parent not exits id=$id, $_";
		next;
	    }
	    
	    my $start = $parent->start;
	    my $csite_start = $start + $spstart*3-2-1  ;
	    my $csite_end = $start + $spend*3-1;  
# need to convert to DNA coordinates
	    my $hasSP = 0;
	    
	    for my $f (@{$seqHash->{$parent->seq_id}{FEATURE}}){
		if($f->primary_tag eq "signal_peptide"){
		    $hasSP = 1;
		    last;
		}
	    }
	    next if $hasSP;
	    $sp_count++;
	    my $sigpep = Bio::SeqFeature::Generic->new(
		-seq_id     => $parent->seq_id,
		-source_tag => "SignalP-4.1",
		-primary    => 'signal_peptide',
		-start      => $csite_start,
		-end        => $csite_end,
		#-score      => $dscore,
		-score      => ".",
		-strand     => $parent->strand, 
		-frame      => '.' # PHASE: compulsory for peptides, can't be '.'
	
		);
	    $sigpep->add_tag_value("Networks",$networks);
	    $sigpep->add_tag_value("D-cutoff", $dcutoff);
	    $sigpep->add_tag_value("ID", "$id\_sp$sp_count");
	    $sigpep->add_tag_value("Parent", $id);
	    $num_sigpep++;
	    push @{$seqHash->{$parent->seq_id}{FEATURE}}, $sigpep; 
	    $parent->add_tag_value("sp","YES");
	}
	
    }
    $/ ="\n";
    
    msg("Found $num_sigpep signal peptides");
    close(SIGNALP);
}


#----------------------------------------------------------------------
## predict transmerbrane helices
sub predict_TM{

    my ($fasta, $seqHash, $seqArray) = @_;

    msg("Looking for transmembrane helices in the predicted proteins");
    my %cds;
        
    for my $sid (@$seqArray) {
	
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    next unless $f->primary_tag eq "CDS";
	    my ($id) = TAG($f, "ID");
	    $cds{$id} = $f;
	    #$cds{++$count} = $f;
	    
	}
    }
    
    my $cmd = "tmhmm --short 1  --workdir $outdir $fasta > $fasta.tmhmm.temp";
    
    if(! -e "$fasta.tmhmm.temp"){
	runcmd($cmd);
    }
    open my $TMHMM, "$fasta.tmhmm.temp";
    
    my $tm_count = 0;
    #open my $TMHMM, '-|', $cmd;
    #METAANNOT_00004 len=509 ExpAA=0.02      First60=0.00    PredHel=0       Topology=o
    #METAANNOT_00005 len=344 ExpAA=213.46    First60=34.76   PredHel=10      Topology=i12-34o49-66i71-88o93-115i124-146o161-183i211-233o243-262i269-291o301-323i
    while (<$TMHMM>) {
	chomp;
	next if /PredHel=0/; 
		
	my @x = split m/\s+/;
	my ($expaa) = $x[2] =~ /ExpAA=(\S+)/;
	my ($first60) = $x[3] =~ /First60=(\S+)/;
	my ($predHel) = $x[4] =~ /PredHel=(\d+)/;
	#convert the coordinate from AA to DNA
	my ($topology) = $x[5] =~ /Topology=(\S+)/;
	my @nums = $topology =~ /(\d+)/g;
	my @chars = $topology =~ /(i|o+)/g;
	my $nc = 0;
	my $i = 0;
	my $cdsid = $x[0];
	my $f = $cds{ $cdsid};
	if(not defined $f){
	    print STDERR "************************************* parent not exits id=$cdsid, $_";
	    next;
	}
	
	my $start = $f->start;
	my $end = $f->end;
	my $strand = $f->strand;
	
	my $dna_topology = "";
	#print STDERR "$_\n";
	#print STDERR join(",", @chars), "\n";
	#print STDERR join(",", @nums), "\n";
	while($i < @chars && $nc < @nums){
	    my $spos = $start + $nums[$nc++]*3-2-1;
	    my $epos = $start + $nums[$nc++]*3-1;
	    $dna_topology .= $chars[$i++] . $spos . "-" . $epos;
	}
	$dna_topology .= "$chars[$#chars]";
#	print STDERR "start=$start , $dna_topology\n";
	$tm_count++;
	my $sid = $f->seq_id;
	my $tmhelix = Bio::SeqFeature::Generic->new(
	    -seq_id     => $sid,
	    -source_tag => "tmhmm",
	    -primary    => 'transmembrane_helix',
	    -start      => $start,
	    -end        => $end,
	    -score      => '.',
	    -strand     => $strand,
	    -frame      => '.'

	    );
	$tmhelix->add_tag_value("ID","$cdsid\_tm\_$tm_count");
	$tmhelix->add_tag_value("Parent",$cdsid);
	$tmhelix->add_tag_value("ExpAA",$expaa);
	$tmhelix->add_tag_value("First60",$first60);
	$tmhelix->add_tag_value("PredHel",$predHel);
	$tmhelix->add_tag_value("Topology", $dna_topology);
	push @{$seqHash->{$sid}{FEATURE}}, $tmhelix; 
	$f->add_tag_value("tm_num",$predHel);

    }
    close($TMHMM);
}

#----------------------------------------------------------------------
# read in sequences
sub init{
    
    my ($contig,$fasta, $mincontiglen, $seqHash, $seqArray) = @_;
    my %seqids = ();
    use FileHandle;
    my $fh = FileHandle->new;
        
    if($contig =~ /.gz$/){
	$fh->open("zcat $contig |")
    }
    else{
	$fh->open($contig)
    }
    my $fout = Bio::SeqIO->new(-file=>">$fasta", -format=>'fasta');
    
    msg( "Loading and checking input file: $contig");
    my $fin = Bio::SeqIO->new(-fh=>$fh, -format=>'fasta');
    my $ncontig = 0;
    while (my $seq = $fin->next_seq) {
	if ($seq->length < $mincontiglen) {
	    msg("Skipping short (<$mincontiglen bp) contig:",$seq->display_id);
	    next;
	}
	$ncontig++;
	my $s = $seq->seq;
	$s = uc($s);
	$s =~ s/[*-]//g;      # replace pads/gaps with nothing
	$s =~ s/[^ACTG]/N/g;  # replace wacky IUPAC with N
	$seq->seq($s);
	$seq->desc(undef);
	
	if (exists $seqids{$seq->id}) {
	    err("Uh oh! Sequence file '$contig' contains duplicate sequence ID:", $seq->id);
	}
	$fout->write_seq($seq);
	$seqids{$seq->id}= 1;
	$seqHash->{ $seq->id }{DNA} = $seq;
	push @$seqArray, $seq->id;  # this array it used to preserve original contig order
    }
    msg("Wrote $ncontig contigs");
    %seqids = ();
}


sub TAG {
  my($f, $tag) = @_;
  return "" unless $f->has_tag($tag);
  my (@values) = ($f->has_tag($tag)) ? $f->get_tag_values($tag) : ("");
  for (@values){
      s/,/%2C/g;
  }
  my $value = join(",", @values);
  
#$value =~ s/^\s+(\S.*?)\s+$/$1/;
  
  return $value;
}

# Option setting routines

sub setOptions {
    use Getopt::Long;
    
    @Options = (
	'General:',
	{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
	{OPT=>"version", VAR=>\&version, DESC=>"Print version and exit"},
	
	'input:',
	{OPT=>"mincontiglen=i",  VAR=>\$mincontiglen, DEFAULT=>200, DESC=>"Minimum contig size [NCBI needs 200]"},
	
	'gene prediction:',
	{OPT=>"gcode=i",  VAR=>\$gcode, DEFAULT=>11, DESC=>"translation table to use for gene predication"},
	{OPT=>"gtype=s",  VAR=>\$gtype, DEFAULT=>"single", DESC=>"single or metagenome: [single, meta]"},
	{OPT=>"minorflen=i",  VAR=>\$minorflen, DEFAULT=>180, DESC=>"Minimum orf length"},
	{OPT=>"evalue=f",VAR=>\$evalue, DEFAULT=>1E-5, DESC=>"evalue cut-off for rRNA prediction"},
	{OPT=>"sp!",  VAR=>\$sp, DEFAULT=>0, DESC=>"Disable signal peptide and cleavage site predication using signalp, it is slow when it is enabled"},
	{OPT=>"tm!",  VAR=>\$tm, DEFAULT=>0, DESC=>"Disable transmembrane helics predication using tmhmm, it is slow when it is enabled"},
	
	'Outputs:',
	{OPT=>"prefix=s",  VAR=>\$prefix, DEFAULT=>'', DESC=>"Filename output prefix"},
	{OPT=>"outdir=s",  VAR=>\$outdir, DEFAULT=>'', DESC=>"Output folder [auto]"},
	{OPT=>"force!",  VAR=>\$force, DEFAULT=>0, DESC=>"Force overwriting existing output folder"},
	{OPT=>"locustag=s",  VAR=>\$locustag, DEFAULT=>$EXE, DESC=>"Locus tag prefix"},	
	'Computation:',
	{OPT=>"cpus=i",  VAR=>\$cpus, DEFAULT=>1, DESC=>"Number of CPUs to use"}
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
	"Name:\n  ", ucfirst($EXE), " $VERSION by $AUTHOR\n",
	"Synopsis:\n  Metagenome contig gene prediction\n",
	"Usage:\n  $EXE [options] <contigs.fasta>\n";
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
sub output_gff{
    
    my ($seqHash, $seqArray, $prefix, $outdir) = @_;
    my $gffver = 3;
    
    msg("Writing all feature gff file to $outdir");
    #open my $gff_fh, '>', "$outdir/$prefix.feature.gff";
    open my $gff_fh, '>', "$outdir/features.gff";
    my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gffver);
    
    print $gff_fh "##gff-version $gffver\n";
    
    for my $id (@seqArray) {
	print $gff_fh "##sequence-region $id 1 ", $seqHash{$id}{DNA}->length, "\n";
    }
    my $idc = 0;
    for my $sid (@seqArray) {
	#my $ctg = $seqHash{$sid}{DNA};
	$idc++;
	for my $f ( sort { $a->start <=> $b->start } @{ $seqHash{$sid}{FEATURE} }) {
	    
	    if($f->has_tag("Parent")){
		my $p = TAG($f, "Parent");
		$p =~ s/$sid/$idc/;
		$f->remove_tag("Parent");
		$f->add_tag_value("Parent", $p);
		
	    }
	    print $gff_fh $f->gff_string($gff_factory),"\n";
	}
    }
    close($gff_fh);
}

sub maskSeqs{
    
    my ($fasta, $output, $locs) = @_;
    
    open(FASTA, $fasta) or die "Could not open $fasta to read, $!\n";
    open(MASKED, ">$output") or die "Could not open $output to write, $!\n";
    
    $/ = "\n>";
    while(<FASTA>){
	chomp;
	if(my ($seqid, $seq) =  /^>?(\S+).*?\n(.*)/s){
	    $seq =~ s/\s+//g;
	    
	    if(exists $locs->{$seqid}){
		foreach my $loc (keys %{$locs->{$seqid}}){
		    my ($start, $end) = split(/:/, $loc);
		    my $flen = $end - $start + 1;
		    substr($seq, $start-1, $flen) = "N"x$flen;
		    #my $flen = $end - $start + 1;
		    #my $front = substr($seq, 0, $start-1);
		    #my $back = substr($seq, $start-1);
		    #$seq = $front . "N"x$flen . $back;
		    
		}
		print MASKED ">$seqid\n$seq\n";
	    }
	    else{
		print MASKED ">$seqid\n$seq\n";
	    }
	}
    }
    $/ = "\n";
    close(FASTA);
    close(MASKED);
}
