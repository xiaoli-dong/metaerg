#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use threads;
use threads::shared;
use Time::Piece;
use Scalar::Util qw(openhandle);
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::SearchIO;
use List::Util qw(min max sum);
#for use fc sort function
use v5.16;
use Benchmark;
use lib "$FindBin::Bin";
require "util.pl";
use List::MoreUtils qw(uniq);

#http://www.genome.jp/kegg-bin/show_pathway?ko00130+K00568+K03182
#download pathway hierachical level: http://www.genome.jp/kegg-bin/get_htext#B1
#http://www.genome.jp/kegg/tool/map_pathway.html
# protein search cutoff: e-5, 30% identity, 70% coverage according JGI metagenome annoation pipeline
#http://www.kegg.jp/kegg/brite.html (Binary Relationships)


#...........................................................................
# Global variabies
my $starttime = localtime;
my $AUTHOR = 'Xiaoli Dong <xdong@ucalgary.ca>';
my $VERSION = "0.1";
my $EXE = $FindBin::RealScript;
my $bin = "$FindBin::RealBin";


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# command line options
my(@Options, $quiet, $debug,$force, $prefix, $outdir,$locustag, $increment, $gcode, $gtype, $mincontiglen, $minorflen, $evalue, $cpus, $identity, $coverage,$hmm_cutoff,$hmm_evalue_cutoff, $sp, $tm, $depth_f, $plevel_f, $DBDIR);

setOptions();
my $sqlite_dir = "$DBDIR/sqlite3";
my $t0 = Benchmark->new;
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# prepare input, output filename, directory
$locustag ||= uc($EXE);
$prefix ||= $locustag.'_'.(localtime->mdy(''));
$outdir ||= $prefix;
my $tmpdir = "$outdir/tmp";
my $datadir = "$outdir/data";

if (-d $outdir) {
   if ($force) {
	msg("Re-using existing --outdir $outdir")
    }
    else {
	err("Folder '$outdir' already exists! Please change --outdir or use --force");
    }
}
else {
    msg("Creating new output folder: $outdir");
    runcmd("mkdir -p \Q$outdir\E");
    runcmd("mkdir -p \Q$tmpdir\E");
    runcmd("mkdir -p \Q$datadir\E");
}


msg("Using filename prefix: $prefix.XXX");

#print out the input parameter values
my $param = "$EXE --dbdir $DBDIR --mincontiglen $mincontiglen --minorflen $minorflen --prefix $prefix --outdir $outdir --locustag $locustag --increment $increment --cpus $cpus --evalue $evalue --identity $identity --coverage $coverage -hmmcutoff $hmm_cutoff -hmmevalue $hmm_evalue_cutoff";

$param .= ($sp) ? " --sp" : "";
$param .= ($tm) ? " --tm" : "";
$param .= ($force) ? " --force" : "";
$param .= " $ARGV[0]\n";
msg($param);


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# prepare fasta, remove short contig, replace ambiguous with Ns
my $contig = shift @ARGV or err("Please supply a contig fasta file on the command line.");
my $fasta = "$outdir/$prefix.fna";
msg("******Start to predicate genes\n");

#contig were filtered in predictFeatures and reformated before gene calling
if(! -e "$fasta"){
    filter_fasta($contig, $fasta, $mincontiglen);
}

my $cmd = "$^X $bin/predictFeatures.pl --dbdir $DBDIR --evalue $evalue --gtype $gtype --gc $gcode --minorflen $minorflen --prefix $prefix --cpus $cpus --outdir $tmpdir ";
$cmd .= "--sp " if ($sp);
$cmd .= "--tm " if ($tm);
$cmd .= "--force $fasta";
runcmd("$cmd");
msg("******Finish predicating genes\n\n");

msg("******Start annotate cds\n");
my $faa = "$tmpdir/cds.faa";
$cmd = "$^X $bin/annotCDs.pl --dbdir $DBDIR --cpus $cpus --evalue $evalue --identity $identity --coverage $coverage --hmmcutoff $hmm_cutoff --hmmevalue $hmm_evalue_cutoff --outdir $tmpdir $faa";
runcmd("$cmd");
msg("******Finish annotate cds\n\n");


my %seqHash = ();

open my $gff, "$tmpdir/features.gff" or die "could not open $tmpdir/features.gff to read, $!\n";

my $gffio = Bio::Tools::GFF->new(-fh =>$gff, -gff_version => 3);

while (my $f = $gffio->next_feature) {
    my $sid = $f->seq_id;
    push (@{$seqHash{$sid}{FEATURE}}, $f);
}
close($gff);

my $fin = Bio::SeqIO->new(-file=>"$fasta", -format=>'fasta');
while (my $seq = $fin->next_seq) {
    $seqHash{$seq->id}{DNA} = $seq;
}

parse_search_results($tmpdir, \%seqHash);

open(MMARKER, "$bin/../txt/metabolic.txt") or die "could not open $bin/../txt/metabolic.txt to read, $!\n";
my %hmm2process = ();

while(<MMARKER>){
    next if /^#/;
    chomp;
    #tr/\r\n//;
    my @l = split(/\t/, $_);
    $hmm2process{$l[3]}->{gene} = $l[2];
    $hmm2process{$l[3]}->{process} = $l[1];
    $hmm2process{$l[3]}->{chemical} = $l[0];
    $hmm2process{$l[3]}->{cutoff_score} = $l[5];

}
close(MMARKER);

use DBI;
my $dbh = DBI->connect(
    "dbi:SQLite:dbname=$sqlite_dir/metaerg.db",
    "",
    "",
    { RaiseError => 1 },
    ) or die $DBI::errstr;

msg("Start mapping sprot db features to cds sequences");
mapping_sprot(\%seqHash);
msg("Finishing mapping sprot db features to cds sequences");


msg("Start mapping metabolic marker gene db features to cds sequences");
mapping_metabolic(\%seqHash);
msg("Finish mapping metabolic marker gene db features to cds sequences");

msg("Start mapping pfam db features to cds sequences");
mapping_pfam(\%seqHash);
msg("Finish mapping pfam db features to cds sequences");

msg("Start mapping tigrfam db features to cds sequences");
mapping_tigrfam(\%seqHash);
msg("Finish mapping tigrfam db features to cds sequences");

msg("Start mapping genomedb db features to cds sequences");
mapping_genomedb(\%seqHash);
msg("Finish mapping genomedb db features to cds sequences");

msg("Start mapping foam db features to cds sequences");
mapping_foam(\%seqHash);
msg("Finish mapping foam db features to cds sequences");


get_all_kos(\%seqHash);
get_all_ecs(\%seqHash);
get_casgene_acc(\%seqHash);
$dbh->disconnect();
#open my $html, ">$outdir/$prefix.html" or die "could not open $outdir/$prefix.html to write, $!\n";
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

###################################################################

# Add locus_tags and protein_id[CDS only] (and Parent genes if asked)
msg("Adding /locus_tag identifiers");
my $num_lt=0;
my %ids2newids = ();
#open (MAP, ">idmap.txt");
for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %seqHash){
	next if not exists $seqHash{$sid}{FEATURE};
	for my $f ( sort { $a->start <=> $b->start } @{ $seqHash{$sid}{FEATURE} }) {
	    next unless $f->primary_tag =~ m/CDS|RNA|signal_peptide|repeat_region|transmembrane_helix/;
	    $num_lt++;
	    my $orgid = ($f->get_tag_values("ID"))[0] if $f->has_tag('ID');
	    my $ID = sprintf("${locustag}|%05d", $num_lt * $increment);
	    $f->add_tag_value('locus_tag', $ID);
	    $f->remove_tag('ID') if $f->has_tag('ID');;
	    $f->add_tag_value('ID', $ID);
	    $ids2newids{$orgid} = $ID;
	    #print MAP "$orgid\t$ID\n";
	}

}
#close(MAP);
msg("Assigned $num_lt locus_tags to CDS and RNA features.");

##########################################################################

msg("Start outputing $prefix gff file");
#if(! -e "$datadir/master.gff"){
    output_gff(\%seqHash,$outdir, \%ids2newids);
#}

msg("Start outputing report files");



$cmd = "$^X $bin/output_reports.pl -g $datadir/master.gff -o $outdir -f $outdir/$prefix.fna";
if($depth_f ne ""){
    $cmd .= " -d $depth_f";
}
if($plevel_f ne ""){
    $cmd .= " -p $plevel_f";
}

runcmd("$cmd");
msg("******Finish outputing report files\n\n");
my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);

msg("metaerg took:" . timestr($td) . " to run\n");

sub parse_search_results{

    my ($outdir, $seqHash) = @_;


    msg("Start parsing diamond search against uniprot_sprot db results");
    my $diamond_output = "$outdir/uniprot_sprot.blasttable";
    parse_diamond($diamond_output,"blasttable","uniprot_sprot", "sprot", \%seqHash, $coverage);
    msg("Finishing parsing diamond search against uniprot_sprot results");


    msg("Start parsing diamond search against genomedb db results");
    $diamond_output = "$outdir/genomedb.blasttable";
    parse_diamond($diamond_output,"blasttable","genomedb", "genomedb", \%seqHash, $coverage);
    msg("Finish parsing diamond search against genomedb db results");

    msg("Start parsing hmmsearch search against pfam db results");
    my $hmmsearch_output = "$outdir/Pfam-A.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","Pfam-A.hmm", "pfam", \%seqHash);
    msg("Finish parsing hmmsearch search against pfam db results");

    msg("Start parsing hmmsearch search against tigrfam db results");
    $hmmsearch_output = "$outdir/TIGRFAMs.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","TIGRFAMs.hmm", "tigrfam", \%seqHash);
    msg("Finish parsing hmmsearch search against tigrfam db results");

    msg("Start parsing hmmsearch search against metabolic_pathway db results");
    $hmmsearch_output = "$outdir/metabolic.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","metabolic.hmm", "homefam", \%seqHash);
    msg("Finish parsing hmmsearch search against metabolic_pathway db results");

    msg("Start parsing hmmsearch search against casgene db results");
    $hmmsearch_output = "$outdir/casgenes.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","casgenes.hmm", "casgene", \%seqHash);
    msg("Finish parsing hmmsearch search against casgenes db results");

    msg("Start parsing hmmsearch search against foam hmm db results");
    $hmmsearch_output = "$outdir/FOAM-hmm_rel1a.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","FOAM-hmm_rel1a.hmm", "foam", \%seqHash);
    msg("Finish parsing hmmsearch search against foam hmm db results");
}

#----------------------------------------------------------------------
sub parse_diamond{

    my ($diamond_output, $format, $dbname, $type, $seqHash, $coverage) = @_;

#qseqid sseqid qlen qstart qend sstart send qframe pident bitscore evalue qcovhsp
    open(F6, $diamond_output) or die "Could not open $diamond_output to read, !$\n";
    my $num_hit = 0;
    while(<F6>){
	chomp;
	my @line = split(/\t/, $_);
	my ($sid) = $line[0] =~ /^(\S+?)\_cds/;
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    if(($f->get_tag_values("ID"))[0] eq $line[0]){
		#my $qlen = ($f->end - $f->start)/3;
		#my $qcov = 100*$a_len/$qlen;
		#if(100*$frac_id >= $identity && $qcov >= $coverage){
		if($line[8] >= $identity && $line[11] >= $coverage){
		    $f->add_tag_value("$type\_target","db:$dbname|$line[1] $line[3] $line[4] evalue:$line[10] qcov:" . sprintf("%.2f", $line[11]). " identity:".sprintf("%.2f", $line[8]));
		    $num_hit++;
		}

	    }
	}
    }
    msg("Found $num_hit $type hits");
    #delfile($diamond_output);
    #return $cds;
}

#----------------------------------------------------------------------
sub parse_diamond_org{

    my ($diamond_output, $format, $dbname, $type, $seqHash) = @_;

    my $bls = Bio::SearchIO->new(-file=>$diamond_output, -format=>$format);
    my $num_hit = 0;
    while (my $res = $bls->next_result) {
	#my $qacc = $res->query_accession;
	#my $qdesc = $res->query_description;
	my $qname = $res->query_name;

	while ( my  $hit = $res->next_hit ){
	    my $hitname = $hit->name;
	    #while ( my $hsp = $hit->next_hsp ){
	    my $hsp = $hit->next_hsp or next;
	    my $sig = $hsp->significance();
	    my $score = $hsp->score;
	    my $start = $hsp->start('query');
	    my $end = $hsp->end('query');
	    my $frac_id = $hsp->frac_identical;
	    my $a_len = $hsp->length("query");
	    #my $f = $cds->{$qname};

	    my ($sid) = $qname =~ /^(\S+?)\_cds/;
	    for my $f (@{$seqHash->{$sid}{FEATURE}}){
		if(($f->get_tag_values("ID"))[0] eq $qname){
		    my $qlen = ($f->end - $f->start)/3;
		    my $qcov = 100*$a_len/$qlen;
		    if(100*$frac_id >= $identity && $qcov >= $coverage){
			$f->add_tag_value("$type\_target","db:$dbname|$hitname $start $end evalue:$sig qcov:" . sprintf("%.2f", $qcov). " identity:".sprintf("%.2f", 100*$frac_id));
			$num_hit++;
		    }

		}
	    }

	}
    }
    msg("Found $num_hit $type hits");
    #delfile($diamond_output);
    #return $cds;
}

sub parse_hmmer3_hmmsearch{
    my ($hmmsearch_out, $format, $dbname, $type, $seqHash) = @_;

    my $hmmer3 = Bio::SearchIO->new(-file=>$hmmsearch_out, -format=>$format);
    #my %dbevalue = ();
    my $num_hit = 0;
    while (my $res = $hmmer3->next_result) {
        my $qacc = $res->query_accession;
        my $qdesc = $res->query_description;
        my $qname = $res->query_name;
        while ( my  $hit = $res->next_hit ){
            my $hitname = $hit->name;
	    my ($sid) = $hitname =~ /^(\S+?)\_cds/;
	    for my $f (@{$seqHash->{$sid}{FEATURE}}){
                if(($f->get_tag_values("ID"))[0] eq $hitname){

		    my $seqT = $hit->score;

		    while ( my $hsp = $hit->next_hsp ){
			my $sig = $hsp->significance();
			my $score = $hsp->score;
			my $start = $hsp->start('hit');
			my $end = $hsp->end('hit');
			my $alignLen = $hsp->length("query");
			my $frac_id = $hsp->frac_identical;
			my $qlen = ($f->end - $f->start)/3;
			my $qcov = 100*$alignLen/$qlen;
			$qacc = $qacc ne "" ? $qacc : $qname;
			my $target_value = "db:$dbname|$qacc $start $end evalue:$sig qcov:" . sprintf("%.2f", $qcov). " identity:".sprintf("%.2f",100*$frac_id) . " score:$score seqT:$seqT name:$qname";
			#print STDERR "qname=$qname, hitname=$hitname, sid=$sid ", "start=", $f->start, " hitname=$hitname\n";
			if($f->has_tag("$type\_target")){
			    my $isOverlap = 0;
			    my $isUpdated = 0;
			    my @values = ();
			    foreach my $ptarget ($f->get_tag_values("$type\_target")){
				my($pstart, $pend, $pevalue, $pscore) = $ptarget =~ /^\S+?\s+(\d+?)\s+(\d+?)\s+?evalue\:(\S+).*?score\:(\S+)/;
				if(is_overlapping($start, $end, $pstart, $pend)){
				    #update the hit which has lower evalue and better score
				    if(($evalue < $pevalue) || ($evalue == $pevalue && $score > $score)){
					push(@values, $target_value);
					$isUpdated = 1;
				    }
				    else{
					push(@values, $ptarget);
				    }

				    $isOverlap = 1;
				}
			    }
			    #not overlap to any existing domains
			    if(!$isOverlap){
				$f->add_tag_value("$type\_target", $target_value);
			    }
			    elsif($isUpdated){
				$f->remove_tag("$type\_target");
				foreach my $value (@values){
				    $f->add_tag_value("$type\_target", $value);
				}
			    }

			}
			else{
			    $f->add_tag_value("$type\_target", $target_value);
			}

		    }
		}
	    }
	}
    }
#delfile($hmmsearch_out);
}
sub get_casgene_acc{
    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {

	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    next unless $f->primary_tag eq "CDS";
	    if($f->has_tag('casgene_target')){

		foreach my $target ($f->get_tag_values('casgene_target')){
		    if(my ($acc) = $target =~ /db:casgenes.hmm\|(\S+)/){

			#for the same gene, multiple domain hit, we only report casgene_acc once for one gene annotation
			if($f->has_tag("casgene_acc")){
			    my $accexist = 0;
			    foreach my $value ($f->get_tag_values("casgene_acc")){

				if($acc eq $value){
				    $accexist = 1;
				    last;
				}
			    }
			    if($accexist == 0){
				$f->add_tag_value("casgene_acc", $acc);
			    }
			}
			else{
			    $f->add_tag_value("casgene_acc", $acc);
			}
		    }
		}
	    }
	}
    }
}
sub get_all_kos{

    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {

	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    #if(($f->get_tag_values("ID"))[0] eq $hitname){

        #for my $fid (keys %{$seqHash->{$sid}}) {
            #product feature
	    #my $f = $seqHash->{$sid}{$fid};
	    #product feature
	    next unless $f->primary_tag eq "CDS";
	    my $foam_kos = ($f->get_tag_values('foam_kos'))[0] if ($f->has_tag('foam_kos'));
	    my $sprot_kos = ($f->get_tag_values('sprot_kos'))[0] if ($f->has_tag('sprot_kos'));

	    $foam_kos //= "";
	    $sprot_kos //= "";
	    my @all_kos = ();

	    push(@all_kos, split(",", $foam_kos)) if $foam_kos ne "";
	    #TODO: update sql database to change the seperate to ,
	    push(@all_kos, split(";", $sprot_kos)) if $sprot_kos ne "";
	    my @all_kos_uniq = uniq(@all_kos) if scalar @all_kos > 0;

	    foreach my $ko (@all_kos_uniq){
		$f->add_tag_value("allko_ids",$ko);
		my $stmt = "SELECT * FROM FOAM_ontology where KO=\"$ko\"";
		my $sth = $dbh->prepare( $stmt);
		$sth->execute();
		my $all = $sth->fetchall_arrayref();
		foreach my $row (@$all) {
		    my ($l1, $l2, $l3, $l4, $ko) = @$row;
		    $f->add_tag_value('allko_ontology', "L1:$l1;L2:$l2;L3:$l3;L4:$l4");
		}
		$sth->finish();
	    }
	}

    }
}

sub get_all_ecs{
    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    #if(($f->get_tag_values("ID"))[0] eq $hitname){

	    #for my $fid (keys %{$seqHash->{$sid}}) {
            #product feature
            #my $f = $seqHash->{$sid}{$fid};
            #product feature
            next unless $f->primary_tag eq "CDS";
            my $foam_ec = ($f->get_tag_values('foam_ec'))[0] if ($f->has_tag('foam_ec'));
            my $sprot_ec = ($f->get_tag_values('sprot_ec'))[0] if ($f->has_tag('sprot_ec'));
	    my $tigrfam_ec = ($f->get_tag_values('tigrfam_ec'))[0] if ($f->has_tag('tigrfam_ec'));
            $foam_ec //= "";
            $sprot_ec //= "";
	    $tigrfam_ec //= "";

            my @all_ecs = ();

            push(@all_ecs, split(",", $foam_ec)) if $foam_ec ne "";
	    push(@all_ecs, split(";", $sprot_ec)) if $sprot_ec ne "";
	    push(@all_ecs, split(";", $tigrfam_ec)) if $tigrfam_ec ne "";
            my @all_ecs_uniq = uniq(@all_ecs) if scalar @all_ecs > 0;

            foreach my $ec (@all_ecs_uniq){
                $f->add_tag_value("allec_ids",$ec);
	    }
        }
    }
}

sub mapping_sprot{
    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    #if(($f->get_tag_values("ID"))[0] eq $hitname){

	    #for my $fid (keys %{$seqHash->{$sid}}) {
            #product feature
	    #my $f = $seqHash->{$sid}{$fid};
	    #product feature
	    next unless $f->primary_tag eq "CDS";
            if ($f->has_tag('sprot_target')) {
		my $target = ($f->get_tag_values('sprot_target'))[0];
		#print "target=$target\n";
		if(my ($spid) = $target =~ /db:uniprot_sprot\|(\S+)/){
		    my $stmt = "select * from uniprot_spid2annot where spid=\"$spid\"";
		    my $sth = $dbh->prepare( $stmt);
		    $sth->execute();
		    my $all = $sth->fetchall_arrayref();

		    foreach my $row (@$all) {
			my ($spid, $ko, $kegg, $ec, $desc,$go, $os, $oc) = @$row;
			#$f->add_tag_value('sprot_OC', $oc) if $oc ne "";
			#$f->add_tag_value('sprot_os', $os) if $os ne "";
			$f->add_tag_value('sprot_desc', $desc) if $desc ne "";
			$f->add_tag_value('sprot_go', $go) if $go ne "";
			$f->add_tag_value('sprot_kos', $ko) if $ko ne "";
			$f->add_tag_value('sprot_kegg', $kegg) if $kegg ne "";
			$f->add_tag_value('sprot_ec', $ec) if $ec ne "";
		    }
		    $sth->finish();


		}
	    }
	}

    }
}

sub mapping_genomedb_org{
    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    #for my $fid (keys %{$seqHash->{$sid}}) {
            #product feature
	    #my $f = $seqHash->{$sid}{$fid};

            #product feature
	    next unless $f->primary_tag eq "CDS";
            if ($f->has_tag('genomedb_target')) {
		#for my $target ($f->get_tag_values('eggNOG_target')){
		my $target = ($f->get_tag_values('genomedb_target'))[0];

		if(my ($taxid, $pid) = $target =~ /db:genome-database\|(\d+?)~~(\S+)/){
		    my $stmt = "select assembly_acc, assembly_level, lineage from genomedb_aa, genomedb_taxid2lineage where geneid=\"$pid\" and genomedb_aa.taxid=genomedb_taxid2lineage.taxid";
		    my $sth = $dbh->prepare( $stmt);
		    $sth->execute();
		    my $all = $sth->fetchall_arrayref();

		    foreach my $row (@$all) {
			my ($assembly_acc, $assembly_level, $lineage) = @$row;
			#$f->add_tag_value('genomedb_assembly_acc',$assembly_acc);
			#$f->add_tag_value('genomedb_assembly_level', $assembly_level);
			$f->add_tag_value('genomedb_OC', $lineage);
			#$f->add_tag_value('genomedb_taxid', $taxid);
		    }
		    $sth->finish();


		}
	    }
	}

    }
}
sub mapping_genomedb{
    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    #for my $fid (keys %{$seqHash->{$sid}}) {
            #product feature
	    #my $f = $seqHash->{$sid}{$fid};

            #product feature
	    next unless $f->primary_tag eq "CDS";
            if ($f->has_tag('genomedb_target')) {
		#for my $target ($f->get_tag_values('eggNOG_target')){
		my $target = ($f->get_tag_values('genomedb_target'))[0];

		if(my ($gid) = $target =~ /db:genomedb\|(\S+?)\|/){
		    my $stmt = "select lineage from genome2taxon where gid=\"$gid\"";
		    
#very slow
#my $stmt = "select lineage from genome2taxon where gid LIKE \"\%$gid\%\"";
		    my $sth = $dbh->prepare( $stmt);
		    $sth->execute();
		    my $all = $sth->fetchall_arrayref();

		    foreach my $row (@$all) {
			my ($lineage) = @$row;
			$lineage =~ s/s__$//;
			#$f->add_tag_value('genomedb_assembly_acc',$assembly_acc);
			#$f->add_tag_value('genomedb_assembly_level', $assembly_level);
			$f->add_tag_value('genomedb_OC', $lineage);
			#$f->add_tag_value('genomedb_taxid', $taxid);
		    }
		    $sth->finish();


		}
	    }
	}

    }
}

sub mapping_pfam{
    my ($seqHash) = @_;
    my $c = 0;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){

        #for my $fid (keys %{$seqHash->{$sid}}) {
            #product feature
	    #my $f = $seqHash->{$sid}{$fid};

	    #product feature
	    next unless $f->primary_tag eq "CDS";


	    if ($f->has_tag('pfam_target')) {

		my %pids = ();
		for my $target ($f->get_tag_values('pfam_target')){
		    if(my ($pfamid, $score) = $target =~ /^db:\S+?\|(\w+).*?\s+.*?seqT:(\S+)/){

			if(not exists $pids{$pfamid}){
			    if(exists $hmm2process{$pfamid}){
				if($hmm2process{$pfamid}->{cutoff_score} eq "-" || $hmm2process{$pfamid}->{cutoff_score} <= $score){
				    #$f->add_tag_value('metabolic_compound',$hmm2process{$pfamid}->{chemical});
				    #$f->add_tag_value('metabolic_process',$hmm2process{$pfamid}->{process});
				    #$f->add_tag_value('metabolic_gene',$hmm2process{$pfamid}->{gene});
				    #$f->add_tag_value('metabolic_geneid',$pfamid);
				    $f->add_tag_value('metabolic_process',"compound:" . $hmm2process{$pfamid}->{chemical} . ";process:" . $hmm2process{$pfamid}->{process} . ";gene:". $hmm2process{$pfamid}->{gene} . ";");
				}
			    }

			    my $stmt = "SELECT acc, id, DE, TP FROM pfams where acc=\"$pfamid\"";
			    my $sth = $dbh->prepare( $stmt);
			    $sth->execute();
			    my $all = $sth->fetchall_arrayref();
			    foreach my $row (@$all) {
				my ($acc, $id, $de, $tp) = @$row;
				$id ||= "UNKNOWN";
				$de ||= "UNKNOWN";
				$tp ||= "UNKNOWN";

				$f->add_tag_value('pfam_desc', $de) if $de ne "";
				$f->add_tag_value('pfam_id', $id) if $id ne "";
				$f->add_tag_value('pfam_type', $tp) if $tp ne "";
			    }
			    $sth->finish();

			    $stmt = "SELECT go_acc from pfam2go where pfam_acc=\"$pfamid\"";
			    $sth = $dbh->prepare( $stmt);
			    $sth->execute();
			    $all = $sth->fetchall_arrayref();
			    my $gos = "";
			    foreach my $row (@$all) {
				my ($go) = @$row;
				$gos .= "$go;";
			    }
			    $gos =~ s/;$//;
			    $f->add_tag_value('pfam_GO', $gos) if $gos ne "";
			    $sth->finish();
			}
			$pids{$pfamid} += 1;;
		    }
		}
	    }

	}
    }

}
sub mapping_tigrfam{
    my ($seqHash) = @_;

    my $c = 0;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){

        #for my $fid (keys %{$seqHash->{$sid}}) {
            #product feature
	    #my $f = $seqHash->{$sid}{$fid};

	    #product feature
	    next unless $f->primary_tag eq "CDS";

	    if ($f->has_tag('tigrfam_target')) {
		for my $target ($f->get_tag_values('tigrfam_target')){
		    if(my ($tigrfamid, $score) = $target =~ /^db:\S+?\|(\S+?)\s+.*?seqT:(\S+)/){
			if(exists $hmm2process{$tigrfamid}){
			    if($hmm2process{$tigrfamid}->{cutoff_score} eq "-" || $hmm2process{$tigrfamid}->{cutoff_score} <= $score){
				$f->add_tag_value('metabolic_process',"compound:" . $hmm2process{$tigrfamid}->{chemical} . ";process:" . $hmm2process{$tigrfamid}->{process} . ";gene:". $hmm2process{$tigrfamid}->{gene} . ";");
			    }
			}
			my $stmt = "SELECT id,de,it,gs,ec,mainrole, sub1role FROM tigrfams where acc=\"$tigrfamid\"";
			my $sth = $dbh->prepare( $stmt);
			$sth->execute();
			my $all = $sth->fetchall_arrayref();
			foreach my $row (@$all) {
			    my ($id,$de,$it,$gs,$ec,$mainrole, $sub1role) = @$row;
			    $de ||= "";
			    $it ||= "";
			    $gs ||= "";
			    $ec ||= "";
			    $mainrole ||= "";
			    $sub1role ||= "";

			    $f->add_tag_value('tigrfam_desc', $de) if $de ne "";
			    $f->add_tag_value('tigrfam_it', $it) if $it ne "";
			    $f->add_tag_value('tigrfam_gs', $gs) if $gs ne "";
			    $f->add_tag_value('tigrfam_ec', $ec) if $ec ne "";
			    $f->add_tag_value('tigrfam_mainrole', $mainrole) if $mainrole ne "";
			    $f->add_tag_value('tigrfam_sub1role', $sub1role) if $sub1role ne "";
			    $f->add_tag_value('tigrfam_id', $id) if $id ne "";

			}
			$sth->finish();

			$stmt = "SELECT go_acc from tigrfam2go where tigrfam_acc=\"$tigrfamid\"";
			$sth = $dbh->prepare( $stmt);
			$sth->execute();
			$all = $sth->fetchall_arrayref();
			my $gos = "";
			foreach my $row (@$all) {
			    my ($go) = @$row;
			    $gos .= "$go;";

			}
			$gos =~ s/;$//;
			$gos ||= "UNKNOWN";
			$f->add_tag_value('tigrfam_GO', $gos) if $gos ne "";
			$sth->finish();

		    }

		}
	    }

	}
    }


}
sub mapping_metabolic{
    my ($seqHash) = @_;

    my $c = 0;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){

        #for my $fid (keys %{$seqHash->{$sid}}) {
            #product feature
	    #my $f = $seqHash->{$sid}{$fid};

	    #product feature
	    next unless $f->primary_tag eq "CDS";
	    if ($f->has_tag('homefam_target')) {
		my %homefams = ();
		for my $target ($f->get_tag_values('homefam_target')){
		    #print STDERR $target, "\n";
		    if(my ($homefamName, $score) = $target =~ /^db:\S+?\|(\S+?)\s+.*?seqT:(\S+)/){
			if(exists $hmm2process{$homefamName}){
			    if(not exists $homefams{$homefamName} && exists $hmm2process{$homefamName}){
				my $cutoff_score = $hmm2process{$homefamName}->{cutoff_score};
				#print STDERR "$homefamName\t$score\tcutoff=$cutoff_score\n";
				if(($cutoff_score eq "-" || $cutoff_score <= $score)
				   #&& $hmm2process{$homefamName}->{chemical} ne ""
				   #&& $hmm2process{$homefamName}->{process} ne ""
				    ){
				    $f->add_tag_value('metabolic_process',"compound:" . $hmm2process{$homefamName}->{chemical} . ";process:" . $hmm2process{$homefamName}->{process} . ";gene:". $hmm2process{$homefamName}->{gene} . ";");

				    $homefams{$homefamName}++;
				}

			    }

			}
			else{
			    $f->remove_tag('homefam_target') if $f->has_tag('homefam_target');
			}

		    }
		}
	    }

	}
    }


}
sub mapping_foam{
    my ($seqHash) = @_;
    my $c = 0;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){

        #for my $fid (keys %{$seqHash->{$sid}}) {
            #product feature
	    #my $f = $seqHash->{$sid}{$fid};
	    #product feature
	    next unless $f->primary_tag eq "CDS";
	    if ($f->has_tag('foam_target')) {

		my %foamids = ();
		my %koids = ();
		my %ecs = ();
		for my $target ($f->get_tag_values('foam_target')){
		    if(my ($foamid,$name) = $target =~ /^\S+?\|(\S+?)\s.*?name:(\S+)/){
			my ($kos, $ecs_str) = split(/\_/, $name);
			$ecs_str //= "";
			$kos //= "";
			if(not exists $foamids{$foamid}){
			    $kos =~ s/KO\://g;
			    map{ $koids{$_}++} split(",", $kos);
			    map{ $ecs{$_}++} split(",", $ecs_str) if $ecs_str ne "";
			}
			$foamids{$foamid} += 1;
		    }
		}

		$f->add_tag_value('foam_kos', join(",", keys %koids));
		$f->add_tag_value('foam_ecs', join(",", keys %ecs)) if keys %ecs > 0;

	    }

	}
    }

}


sub output_gff{

    my ($seqHash,$outdir, $ids2newids) = @_;
    my $gffver = 3;
    msg("Writing master.gff file to $outdir/data");
    open my $gff_fh, '>', "$outdir/data/master.gff";
    #print STDERR join("***********\n", keys %$ids2newids);
    my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gffver);

    print $gff_fh "##gff-version $gffver\n";

    for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %$seqHash) {
	for my $f ( sort { $a->start <=> $b->start } @{ $seqHash->{$sid}{FEATURE} }) {
	    print $gff_fh $f->gff_string($gff_factory),"\n";
        }
    }
    print $gff_fh "##FASTA\n";
    my $fin = Bio::SeqIO->new(-file=>"$outdir/tmp/cds.faa", -format=>'fasta');
    while (my $seq = $fin->next_seq) {

	print $gff_fh ">", $ids2newids->{$seq->id}, "\n", $seq->seq(), "\n";

    }
    close($gff_fh);
}


# Option setting routines

sub setOptions {
    use Getopt::Long;

    @Options = (
	'General:',
	{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
	{OPT=>"version", VAR=>\&version, DESC=>"Print version and exit"},

	'input:',
	{OPT=>"dbdir=s",  VAR=>\$DBDIR, DEFAULT=>"db", DESC=>"metaerg searching database directory"},
	{OPT=>"mincontiglen=i",  VAR=>\$mincontiglen, DEFAULT=>200, DESC=>"Minimum contig size [NCBI needs 200]"},

	'gene prediction:',
	{OPT=>"gcode=i",  VAR=>\$gcode, DEFAULT=>11, DESC=>"translation table to use for gene predication"},
	{OPT=>"gtype=s",  VAR=>\$gtype, DEFAULT=>"meta", DESC=>"single or metagenome: [arc|bac|euk|meta]"},
	{OPT=>"minorflen=i",  VAR=>\$minorflen, DEFAULT=>180, DESC=>"Minimum orf length"},
	{OPT=>"sp!",  VAR=>\$sp, DEFAULT=>0, DESC=>"Enable signal peptide and cleavage site predication using signalp, it is slow when it is enabled"},
	{OPT=>"tm!",  VAR=>\$tm, DEFAULT=>0, DESC=>"Enable transmembrane helics predication using tmhmm, it is slow when it is enabled"},

	'Outputs:',
	{OPT=>"prefix=s",  VAR=>\$prefix, DEFAULT=>'', DESC=>"Filename output prefix"},
	{OPT=>"outdir=s",  VAR=>\$outdir, DEFAULT=>'', DESC=>"Output folder [auto]"},
	{OPT=>"force!",  VAR=>\$force, DEFAULT=>0, DESC=>"Force overwriting existing output folder"},
	{OPT=>"locustag=s",  VAR=>\$locustag, DEFAULT=>$EXE, DESC=>"Locus tag prefix"},
	{OPT=>"increment=i",  VAR=>\$increment, DEFAULT=>1, DESC=>"Locus tag counter increment"},

	'Computation:',
	{OPT=>"cpus=i",  VAR=>\$cpus, DEFAULT=>8, DESC=>"Number of CPUs to use"},

	'diamond cutoff:',
	{OPT=>"evalue=f",  VAR=>\$evalue, DEFAULT=>1e-5, DESC=>"Similarity e-value cut-off"},
	{OPT=>"identity=f",  VAR=>\$identity, DEFAULT=>20, DESC=>"identity"},
	{OPT=>"coverage=f",  VAR=>\$coverage, DEFAULT=>70, DESC=>"coverage"},

	'hmmsearch cutoff:',
	{OPT=>"hmmcutoff=s",  VAR=>\$hmm_cutoff, DEFAULT=>"--cut_tc", DESC=>"hmm search trusted score threshold: [--cut_ga|--cut_nc|--cut_tc]"},
	{OPT=>"hmmevalue=f",  VAR=>\$hmm_evalue_cutoff, DEFAULT=>"1e-5", DESC=>"custom hmm db search threshold"},

	'abundance input:',
	{OPT=>"depth=s",  VAR=>\$depth_f, DEFAULT=>"", DESC=>"contig coverage depth file in format of:contigid contigLe depth"},
	{OPT=>"plevel=s",  VAR=>\$plevel_f, DEFAULT=>"", DESC=>"protein expession level input file in format of: cds_id cds_desc expression_level"}
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
# read in sequences; remove small contigs; replace ambig with N
sub filter_fasta{

    my ($contig, $fasta, $mincontiglen) = @_;
    my %seqids = ();
    msg( "Loading and checking input file: $fasta");
    my $fin = Bio::SeqIO->new(-file=>"$contig", -format=>'fasta');
    my $fout = Bio::SeqIO->new(-file=>">$fasta", -format=>'fasta');
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

    }
    $ncontig > 0 or err("FASTA file '$contig' contains no suitable sequence entries");
    msg("Wrote $ncontig contigs");
    %seqids = ();
}


#----------------------------------------------------------------------
sub usage {
    print STDERR
	"Name:\n  ", ucfirst($EXE), " $VERSION by $AUTHOR\n",
	"\nSynopsis:\n  Metagenome annotation\n",
	"\nUsage:\n  $EXE [options] <contigs.fasta>\n";
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
	    print STDERR "\n$_\n";
	}
    }
    exit(1);
 }

#----------------------------------------------------------------------

sub version {
    print STDERR "$EXE $VERSION\n";
    exit;
}

