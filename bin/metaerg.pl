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
my $VERSION = "1.2.0";
my $EXE = $FindBin::RealScript;
my $bin = "$FindBin::RealBin";

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# command line options
my(@Options, $quiet, $debug,$force, $prefix, $outdir,$locustag, $increment, $gcode, $gtype, $mincontiglen, $minorflen, $evalue, $cpus, $identity, $coverage,$hmm_cutoff,$hmm_evalue_cutoff, $sp, $tm, $depth_f, $DBDIR);

setOptions();
$DBDIR ||= "$bin/../db";
my $sqlite_dir = "$DBDIR/sqlite3";
my $t0 = Benchmark->new;
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# prepare input, output filename, directory
$locustag ||= uc($EXE);
$prefix ||= $locustag.'_'.(localtime->mdy(''));
$outdir ||= $prefix;
my $tmpdir = "$outdir/tmp";
my $datadir = "$outdir/data";
my $txtdir = "$bin/../txt";
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

my %contig2depth = ();
if($depth_f ne ""){
    my $ref = parse_depth_file($depth_f);
    %contig2depth = %$ref;
}


#contig were filtered in predictFeatures and reformated before gene calling
if(! -e "$fasta"){
    filter_fasta($contig, $fasta, $mincontiglen);
}
msg("******Start to predicate genes\n");
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

open my $gff, "$tmpdir/features.annot.gff" or die "could not open $tmpdir/features.annot.gff to read, $!\n";

my $gffio = Bio::Tools::GFF->new(-fh =>$gff, -gff_version => 3);

while (my $f = $gffio->next_feature) {
    my $sid = $f->seq_id;

    if(exists $contig2depth{$sid}){
	my $depth_names = ();
	my $depth_values = ();
	for my $key (keys %{$contig2depth{$sid}}){
	    $f->add_tag_value("mdepth_cols", $key);
	    $f->add_tag_value("mdepth_values", $contig2depth{$sid}->{$key});
	}
    }
    push (@{$seqHash{$sid}{FEATURE}}, $f);
}
close($gff);

my $fin = Bio::SeqIO->new(-file=>"$fasta", -format=>'fasta');
while (my $seq = $fin->next_seq) {
    $seqHash{$seq->id}{DNA} = $seq;
}

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
	    #$f->add_tag_value('locus_tag', $ID);
	    $f->remove_tag('ID') if $f->has_tag('ID');
	    $f->add_tag_value('ID', $ID);
	    $ids2newids{$orgid} = $ID;
	    #print MAP "$orgid\t$ID\n";
	}

}
#close(MAP);
msg("Assigned $num_lt locus_tags to CDS and RNA features.");

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

get_all_gos(\%seqHash);
get_casgene_acc(\%seqHash);

$dbh->disconnect();


##########################################################################

msg("Start outputing $prefix gff file");
#if(! -e "$datadir/master.gff"){
    output_gff(\%seqHash,$outdir, \%ids2newids);
#}

msg("Start outputing report files");



$cmd = "$^X $bin/output_reports.pl -g $datadir/all.gff -o $outdir -f $outdir/$prefix.fna -db $DBDIR";


runcmd("$cmd");
msg("******Finish outputing report files\n\n");

msg("Start creating result package");
$cmd = "tar --exclude=\'$outdir/tmp\' -cvzf $outdir.tar.gz $outdir";
runcmd("$cmd");
msg("******Finish creating the result package\n\n");

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);

msg("metaerg took:" . timestr($td) . " to run\n");

sub get_all_kos{

    my ($seqHash) = @_;

    for my $sid (keys %$seqHash) {

	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    my $featureid = ($f->get_tag_values("ID"))[0];
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
	    }
	}

    }


}


sub get_all_ecs{
    my ($seqHash) = @_;

    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    my $featureid = ($f->get_tag_values("ID"))[0];
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
	    push(@all_ecs, split(/\s+/, $tigrfam_ec)) if $tigrfam_ec ne "";
            my @all_ecs_uniq = uniq(@all_ecs) if scalar @all_ecs > 0;

            foreach my $ec (@all_ecs_uniq){
                $f->add_tag_value("allec_ids",$ec);

	    }
        }
    }

}


sub get_all_gos{
    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    next unless $f->primary_tag eq "CDS";
            my $pfam_go = $f->has_tag('pfam_go') ? ($f->get_tag_values('pfam_go'))[0] : "";
            my $sprot_go = $f->has_tag('sprot_go') ? ($f->get_tag_values('sprot_go'))[0] : "";
            my @all_gos = ();

            push(@all_gos, split(";", $pfam_go)) if $pfam_go ne "";
	    push(@all_gos, split(";", $sprot_go)) if $sprot_go ne "";
	    my @all_gos_uniq = uniq(@all_gos) if scalar @all_gos > 0;

            foreach my $go (@all_gos_uniq){
                $f->add_tag_value("allgo_ids",$go);
	    }
        }
    }
}

sub mapping_sprot{
    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){

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
			$f->add_tag_value('sprot_id', $spid) if $spid ne "";
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


sub mapping_genomedb{
    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){

	    next unless $f->primary_tag eq "CDS";
            if ($f->has_tag('genomedb_target')) {
		#for my $target ($f->get_tag_values('eggNOG_target')){
		my $target = ($f->get_tag_values('genomedb_target'))[0];

		if(my ($gid) = $target =~ /db:genomedb\|(\S+?)\|/){
		    $f->add_tag_value('genomedb_acc', $gid);
		    my $stmt = "select lineage from genome2taxon where gid=\"$gid\"";

#very slow
#my $stmt = "select lineage from genome2taxon where gid LIKE \"\%$gid\%\"";
		    my $sth = $dbh->prepare( $stmt);
		    $sth->execute();
		    my $all = $sth->fetchall_arrayref();

		    foreach my $row (@$all) {
			my ($lineage) = @$row;
			$lineage =~ s/s__$//;
			$f->add_tag_value('genomedb_OC', $lineage);

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

	    next unless $f->primary_tag eq "CDS";
	    if ($f->has_tag('pfam_target')) {

		for my $target ($f->get_tag_values('pfam_target')){
		    if(my ($pfamid) = $target =~ /^db:\S+?\|(\w+).*/){
			my $stmt = "SELECT acc, id, DE, TP FROM pfams where acc=\"$pfamid\"";
			my $sth = $dbh->prepare( $stmt);
			$sth->execute();
			my $all = $sth->fetchall_arrayref();
			foreach my $row (@$all) {
			    my ($acc, $id, $de, $tp) = @$row;
			    $id ||= "UNKNOWN";
			    $de ||= "UNKNOWN";
			    $tp ||= "UNKNOWN";
			    $f->add_tag_value('pfam_acc', $acc) if $acc ne "";
			    $f->add_tag_value('pfam_desc', $de) if $de ne "";
			    $f->add_tag_value('pfam_id', $id) if $id ne "";
			    #$f->add_tag_value('pfam_type', $tp) if $tp ne "";
			    #print STDERR "pfam=$target\n";
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
			$f->add_tag_value('pfam_go', $gos) if $gos ne "";
			$sth->finish();


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

	    next unless $f->primary_tag eq "CDS";

	    if ($f->has_tag('tigrfam_target')) {
		for my $target ($f->get_tag_values('tigrfam_target')){
		    if(my ($tigrfamid) = $target =~ /^db:\S+?\|(\S+)/){
			my $stmt = "SELECT acc,id,de,it,gs,ec,mainrole, sub1role FROM tigrfams where acc=\"$tigrfamid\"";
			my $sth = $dbh->prepare( $stmt);
			$sth->execute();
			my $all = $sth->fetchall_arrayref();
			foreach my $row (@$all) {
			    my ($acc, $id,$de,$it,$gs,$ec,$mainrole, $sub1role) = @$row;
			    $de ||= "";
			    $it ||= "";
			    $gs ||= "";
			    $ec ||= "";
			    $mainrole ||= "";
			    $sub1role ||= "";
			    $f->add_tag_value('tigrfam_acc', $acc) if $acc ne "";
			    $f->add_tag_value('tigrfam_desc', $de) if $de ne "";
			    #$f->add_tag_value('tigrfam_it', $it) if $it ne "";
			    #$f->add_tag_value('tigrfam_gs', $gs) if $gs ne "";
			    $f->add_tag_value('tigrfam_ec', $ec) if $ec ne "";
			    $f->add_tag_value('tigrfam_mainrole', $mainrole) if $mainrole ne "";
			    $f->add_tag_value('tigrfam_sub1role', $sub1role) if $sub1role ne "";
			    $f->add_tag_value('tigrfam_name', $id) if $id ne "";

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
			#$gos ||= "UNKNOWN";
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

    my $c = 0;
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){

	    next unless $f->primary_tag eq "CDS";
	    if ($f->has_tag('metabolic_target')) {
		my %metabolicfams = ();
		#print STDERR "xiaoli********************************************\n";
		for my $target ($f->get_tag_values('metabolic_target')){
		    #print STDERR $target, "\n";
		    if(my ($metabolicfamName, $score) = $target =~ /^db:\S+?\|(\S+?)\s+.*?score:(\S+)/){
			#print STDERR "xiaoli if metabolicfamName=$metabolicfamName, score=$score********************************************\n";
			if(exists $hmm2process{$metabolicfamName}){
			    if(not exists $metabolicfams{$metabolicfamName} && exists $hmm2process{$metabolicfamName}){
				my $cutoff_score = $hmm2process{$metabolicfamName}->{cutoff_score};
				print STDERR "$metabolicfamName\t$score\tcutoff=$cutoff_score\n";
				if(($cutoff_score eq "-" || $cutoff_score <= $score)){
				    $f->add_tag_value('metabolic_process',"compound:" . $hmm2process{$metabolicfamName}->{chemical} . ";process:" . $hmm2process{$metabolicfamName}->{process} . ";gene:". $hmm2process{$metabolicfamName}->{gene} . ";");
				    #$f->add_tag_value('metabolic_process',"compound:" . $hmm2process{$metabolicfamName}->{chemical} . ";process:" . $hmm2process{$metabolicfamName}->{process} .  ";");
				    $f->add_tag_value("metabolic_acc",$metabolicfamName);
				    $metabolicfams{$metabolicfamName}++;
				}
				else{
				    #print STDERR "************************1 remove $metabolicfamName\n";
				    $f->remove_tag('metabolic_target') if $f->has_tag('metabolic_target');
				}


			    }
			    else{
				#print STDERR "************************2remove $metabolicfamName\n";
				$f->remove_tag('metabolic_target') if $f->has_tag('metabolic_target');
			    }

			}
			else{
			    #print STDERR "************************3remove $metabolicfamName\n";
			    $f->remove_tag('metabolic_target') if $f->has_tag('metabolic_target');
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

	    next unless $f->primary_tag eq "CDS";
	    if ($f->has_tag('foam_target')) {

		#my %foamids = ();
		my %koids = ();
		my %ecs = ();
		for my $target ($f->get_tag_values('foam_target')){
		    if(my ($acc,$name) = $target =~ /^\S+?\|(\S+?)\s.*?name:(\S+)/){

			#$f->add_tag_value('foam_acc', $acc) if $acc ne "";
			my ($kos, $ecs_str) = split(/\_/, $name);
			$ecs_str //= "";
			$kos //= "";
			#if(not exists $foamids{$acc}){
			$kos =~ s/KO\://g;
			map{ $koids{$_}++} split(",", $kos);
			map{ $ecs{$_}++} split(",", $ecs_str) if $ecs_str ne "";
			#}
			#$foamids{$acc} += 1;
		    }
		}

		$f->add_tag_value('foam_kos', join(",", keys %koids));
		$f->add_tag_value('foam_ecs', join(",", keys %ecs)) if keys %ecs > 0;

	    }

	}
    }

}
sub get_casgene_acc{
    my ($seqHash) = @_;
    for my $sid (keys %$seqHash) {

	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    next unless $f->primary_tag eq "CDS";
	    if($f->has_tag('casgene_target')){

		foreach my $target ($f->get_tag_values('casgene_target')){
		    if(my ($acc, $full_score, $best_domain_score, $name) = $target =~ /db:casgenes.hmm\|(\S+?)\s.*?score:(\S+).*?best_domain_score:(\S+?)\s+name:(\S+)/){


			#for the same gene, multiple domain hit, we only report casgene_acc once for one gene annotation
			$f->add_tag_value("casgene_acc", $acc);
			$f->add_tag_value("casgene_name", $name);

		    }
		}

	    }
	}
    }
}


sub output_gff{

    my ($seqHash,$outdir, $ids2newids) = @_;
    my $gffver = 3;
    msg("Writing all.gff file to $outdir/data");
    open my $all_gff_fh, '>', "$outdir/data/all.gff";

   # msg("Writing master.gff.txt file to $outdir/data");
    #open my $master_gff_fh, '>', "$outdir/data/master.gff.txt";

    my $all_gff_factory = Bio::Tools::GFF->new(-gff_version=>$gffver);
    #my $master_gff_factory = Bio::Tools::GFF->new(-gff_version=>$gffver);
    print $all_gff_fh "##gff-version $gffver\n";
    #print $master_gff_fh "##gff-version $gffver\n";

    for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %$seqHash) {
	for my $f ( sort { $a->start <=> $b->start } @{ $seqHash->{$sid}{FEATURE} }) {

	    if($f->primary_tag =~ m/signal_peptide|transmembrane_helix/){
		$f->remove_tag("Parent");

	    }
	    print $all_gff_fh $f->gff_string($all_gff_factory),"\n";

	    foreach my  $ftag ($f->all_tags()){
		if($ftag =~ /foam_ecs|foam_kos|foam_target|sprot_kegg|sprot_kos|sprot_ec|sprot_go|pfam_go|tigrfam_ec|tigrfam_mainrole|tigrfam_sub1role|_target/){
		    $f->remove_tag($ftag);
		}
	    }
	    #print $master_gff_fh $f->gff_string($master_gff_factory),"\n";
        }
    }
    print $all_gff_fh "##FASTA\n";
    my $fin = Bio::SeqIO->new(-file=>"$outdir/tmp/cds.faa", -format=>'fasta');
    while (my $seq = $fin->next_seq) {
	print $all_gff_fh ">", $ids2newids->{$seq->id}, "\n", $seq->seq(), "\n";

    }
    close($all_gff_fh);
    #close($master_gff_fh);
}
sub parse_depth_file{

    my ($fdepth) = @_;
    my %contig2depth = ();
    my %headerIndex2SampleName = ();

    open FDEPTH, $fdepth or die "could not open $fdepth to read, $!\n";
    #header: contigName	contigLen	totalAvgDepth	S1.bam	S1.bam-var	S2.bam	S2.bam-var	S3.bam	S3.bam-var
    while (<FDEPTH>) {
	chomp;
	next if /^#/;
	#header containing sample name info
	if(/^contigName/){
	    my @header = split(/\t/, $_);
	    for (my $i = 3; $i < @header; $i++){
		next if $header[$i] =~ /-var$/;
		my $sampleName = $header[$i];
		$sampleName =~ s/\.bam|_sorted//g;
		$headerIndex2SampleName{$i} = $sampleName;
	    }
	}
	else{
	    my @l = split(/\t/, $_);
	    my ($contigName) = $l[0] =~ /^(\S+)/;
	    my $totalAvgDepth = $l[2];

            #read in Depth
	    $contig2depth{$contigName}->{totalAvgDepth} = $totalAvgDepth;
	    for my $i (keys %headerIndex2SampleName){
		$contig2depth{$contigName}->{$headerIndex2SampleName{$i}} = $l[$i];
	    }
	}
    }
    close(FDEPTH);

    return \%contig2depth;
}


# Option setting routines

sub setOptions {
    use Getopt::Long;

    @Options = (
	'General:',
	{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
	{OPT=>"version", VAR=>\&version, DESC=>"Print version and exit"},

	'input:',
	{OPT=>"dbdir=s",  VAR=>\$DBDIR, DEFAULT=>'', DESC=>"metaerg searching database directory"},
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
	{OPT=>"hmmcutoff=s",  VAR=>\$hmm_cutoff, DEFAULT=>"--cut_ga", DESC=>"hmm search trusted score threshold: [--cut_ga|--cut_nc|--cut_tc]"},
	{OPT=>"hmmevalue=f",  VAR=>\$hmm_evalue_cutoff, DEFAULT=>"1e-5", DESC=>"custom hmm db search threshold"},

	'abundance input:',
	{OPT=>"depth=s",  VAR=>\$depth_f, DEFAULT=>"", DESC=>"contig coverage depth file in format of:contigid contigLe depth"}
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

