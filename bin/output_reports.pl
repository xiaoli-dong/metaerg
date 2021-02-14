#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use HTML::Entities;
use Bio::Tools::GFF;
use Bio::SeqIO;
use List::Util qw(min max sum);
use FindBin;
use lib "$FindBin::Bin";
require "util.pl";
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use Bio::SeqFeature::Generic;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::SimpleValue;
use Data::Dumper;
use File::Path qw(make_path remove_tree);


my ($gff,$fasta,$DBDIR);
my $outdir = ".";
my $bin = "$FindBin::RealBin";
my $templatedir = "$bin/../template";
my $EXE = $FindBin::RealScript;
my $txtdir = "$bin/../txt";
&GetOptions(
    "g=s" =>\$gff,
    "f=s" =>\$fasta,
    "o=s" =>\$outdir,
	"db=s" => \$DBDIR
    );
my $sqlite3_dir = "$DBDIR/sqlite3";

($gff && $fasta && $outdir) ||
    die "Name:\n".
    "$EXE by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "  Generate all all type of reports based on the features in the gff file\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -g <input gff file: all.gff file produced by metaerg>\n".
    "  -f <fasta format contig file, the features were assocaited with this input file>\n".
    "  -o <output dir>\n";
    
my $datadir = "$outdir/data";

if (-d $datadir) {
    msg("Re-using existing --outdir $outdir")
}
else {
    msg("Creating new output folder: $outdir");
    runcmd("mkdir -p \Q$outdir\E");
    runcmd("mkdir -p \Q$datadir\E");
}
use DBI;
my $dbh = DBI->connect(
    "dbi:SQLite:dbname=$sqlite3_dir/metaerg.db",
    "",
    "",
    { RaiseError => 1 },
    ) or die $DBI::errstr;



#print STDERR join("proteomics *************\n", keys %cds2plevels);
my %seqHash = ();

#read in feature.gff file
msg("start to read $gff file");
open my $gff_handle, $gff or die "could not open $gff to read, $!\n";
my $gffio = Bio::Tools::GFF->new(-fh =>$gff_handle, -gff_version => 3);


while (my $f = $gffio->next_feature) {
    my $sid = $f->seq_id;
    push (@{$seqHash{$sid}{FEATURE}}, $f);
}
my %gene2pathways = ();
my %cds_aa_seqs = map { $_->id => $_ } $gffio->get_seqs();
close($gff_handle);
msg("end to read $gff file");
msg("start to read $fasta file");

my $fin = Bio::SeqIO->new(-file=>"$fasta", -format=>'fasta');
while (my $seq = $fin->next_seq) {
    $seqHash{$seq->id}{DNA} = $seq;
}
msg("end to read $fasta file");

my ($ko2genes, $ec2genes) = output_geneAnnotation(\%seqHash, $datadir);
my $keggs = predict_kegg_pathways("$datadir/cds.gene2ko", $ko2genes, \%seqHash);
my $metacycs = predict_metacyc_pathways("$datadir/cds.gene2ec", $ec2genes, \%seqHash);

output_master_annot_summary(\%seqHash, $datadir);
output_tbl(\%seqHash, $datadir);
output_stats(\%seqHash, $datadir);
output_fasta(\%seqHash, \%cds_aa_seqs, $datadir);
output_profiles(\%seqHash, $datadir);
output_gff(\%seqHash,$datadir);
output_htmlreport($templatedir, "$outdir");
$dbh->disconnect();
sub output_htmlreport{

    my($templatedir, $outdir) = @_;

    remove_tree "$outdir/html" if -d "$outdir/html";
    dircopy("$templatedir/html","$outdir/html", ) or die "Copy $templatedir/html to $outdir/html failed, $!\n";
    dircopy("$templatedir/images","$outdir/images") or die "Copy $templatedir/images to $outdir/images failed, $!\n";
    dircopy("$templatedir/js","$outdir/js") or die "Copy $templatedir/js to $outdir/js failed, $!\n";
    copy "$templatedir/help.html", "$outdir"  or die "Copy $templatedir/help.html to $outdir/help.html failed, $!\n";
    copy "$templatedir/index.html", "$outdir"  or die "Copy $templatedir/index.html to $outdir/index.html failed, $!\n";
    copy "$templatedir/style.css", "$outdir" or die "Copy $templatedir/style.css to $outdir/style.css failed, $!\n";
   
    
}

sub output_master_annot_summary{
    my ($seqHash,$datadir) = @_;
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# Write it all out!
    msg("Writing tabular format summary output to $datadir/");
    open my $tbl_fh, '>', "$datadir/master.tsv.txt";
    my @tags = ("casgene_acc", "sp", "tm_num", "sprot_desc", "tigrfam_desc", "pfam_desc","genomedb_OC");

    print $tbl_fh "#contigid\t";
    print $tbl_fh "feature_id\t";
    print $tbl_fh "type\t";
    print $tbl_fh "start\t";
    print $tbl_fh "end\t";
    print $tbl_fh "length\t";
    print $tbl_fh "strand\t";
    print $tbl_fh "kegg_pathways\t";
    print $tbl_fh "metacyc_pathways\t";
    print $tbl_fh join("\t", @tags), "\n";

    for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %$seqHash) {
	next if not exists $seqHash->{$sid}{FEATURE};
        for my $f ( sort { $a->start <=> $b->start } @{ $seqHash->{$sid}{FEATURE} }) {
	    next if $f->primary_tag eq "transmembrane_helix";
            next if $f->primary_tag eq "signal_peptide";
	    my $geneid = ($f->get_tag_values("ID"))[0] if $f->has_tag("ID");
	    my @kegg_pathways = ();
	    my @metacyc_pathways = ();
	    if(exists $gene2pathways{$geneid}->{KEGG}){
		push (@kegg_pathways, keys %{$gene2pathways{$geneid}->{KEGG}});
	    }
	    if(exists $gene2pathways{$geneid}->{metacyc}){
		push (@metacyc_pathways, keys %{$gene2pathways{$geneid}->{metacyc}});
	    }
	    print $tbl_fh "$sid\t";
            print $tbl_fh ($f->get_tag_values("ID"))[0], "\t" if $f->has_tag("ID");
            print $tbl_fh $f->primary_tag, "\t";
            print $tbl_fh $f->start, "\t";
            print $tbl_fh $f->end, "\t";
	    print $tbl_fh $f->length, "\t";
            print $tbl_fh $f->strand, "\t";
	    print $tbl_fh join(";", map{"ko$_"} @kegg_pathways), "\t";
	    print $tbl_fh join(";", @metacyc_pathways), "\t";
            foreach my $tag (@tags){
                print $tbl_fh TAG($f, $tag), "\t";
            }
            print $tbl_fh "\n";
        }
    }

    close($tbl_fh);
}

sub output_stats{
    my ($seqHash, $datadir) = @_;
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# Output a general .txt file with statistics about the annotation

    msg("Generating annotation statistics file");
    my @contig_lens = map {$seqHash{$_}{DNA}->length } keys %$seqHash;
    my @sorted_lens = sort { $a <=> $b} @contig_lens;

    open my $stat_fh, '>', "$datadir/master.stats.txt";
    print $stat_fh "##############input contig stats#############\n\n";
    print $stat_fh "Contig count: ", scalar(@contig_lens), "\n";
    print $stat_fh "Contig total bases: ", sum(@contig_lens), "\n";
    print $stat_fh "Contig min length: ", $sorted_lens[0], "\n";
    print $stat_fh "Contig max length: ", $sorted_lens[$#sorted_lens], "\n";
    print $stat_fh "Contig mean length: ", mean(\@contig_lens), "\n";
    print $stat_fh "Contig meadian length: ", median(\@contig_lens), "\n";
    print $stat_fh "Contig length stdev: ", stdev(\@contig_lens), "\n";
    print $stat_fh "Contig N50: ", get_N50(\@contig_lens), "\n";

    

    my %count;
    for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %$seqHash) {
	next if not exists $seqHash->{$sid}{FEATURE};
        for my $f (@{ $seqHash{$sid}{FEATURE} }) {
	    $count{ $f->primary_tag}->{count}++;
	    $count{$f->primary_tag}->{foam}++ if $f->has_tag('foam_kos');
            $count{$f->primary_tag}->{tigrfam}++ if $f->has_tag('tigrfam_desc');
            $count{$f->primary_tag}->{pfam}++ if $f->has_tag('pfam_desc');
            $count{$f->primary_tag}->{sprot}++ if $f->has_tag('sprot_desc');
            $count{$f->primary_tag}->{genomedb}++ if $f->has_tag('genomedb_OC');
            $count{$f->primary_tag}->{metabolic}++ if $f->has_tag('metabolic_process');
	    $count{$f->primary_tag}->{rRNA_taxon}++ if $f->has_tag('rRNA_taxon');
        }
    }

    my $gene_predict_str = "";
    my $gene_annot_str = "";
    
    for my $ft (sort keys %count) {
	#print $stat_fh "$ft count:\t", $count{$ft}->{count}, "\n";
	$gene_predict_str .= "$ft count:\t". $count{$ft}->{count}. "\n";
	
	if($ft eq "CDS"){
	    $gene_annot_str .= (not exists $count{$ft}->{sprot}) ? "$ft swissprot_hit:\t0\n"  :  "$ft swissprot_hit:\t" . $count{$ft}->{sprot} . "\n";
	    $gene_annot_str .= (not exists $count{$ft}->{pfam}) ? "$ft pfam_hit:\t0\n"  :  "$ft pfam_hit:\t" . $count{$ft}->{pfam} . "\n";
	    $gene_annot_str .= (not exists $count{$ft}->{tigrfam}) ? "$ft tigrfam_hit:\t0\n"  :  "$ft tigrfam_hit:\t" . $count{$ft}->{tigrfam} . "\n";
	    $gene_annot_str .= (not exists $count{$ft}->{foam}) ? "$ft foam_hit:\t0\n"  :  "$ft foam_hit:\t" . $count{$ft}->{foam} . "\n";
	    $gene_annot_str .= (not exists $count{$ft}->{metabolic}) ? "$ft metabolic_marker_hit:\t0\n"  :  "$ft metabolic_marker_hit:\t" . $count{$ft}->{metabolic} . "\n";
	    $gene_annot_str .= (not exists $count{$ft}->{genomedb}) ? "$ft genomedb_hit:\t0\n"  :  "$ft genomedb_hit:\t" . $count{$ft}->{genomedb} . "\n";
	}
	elsif($ft =~ /SrRNA/){
	    $gene_annot_str .= (exists $count{$ft}->{rRNA_taxon}) ?  "$ft SILVA database search hits:\t" . $count{$ft}->{rRNA_taxon} . "\n" : "$ft silva database searching hits:\t0\n";
	}
	
    }
    print $stat_fh "##############gene prediction stats#############\n";
    print $stat_fh $gene_predict_str;
    print $stat_fh "##############gene annotation stats#############\n";
    print $stat_fh $gene_annot_str;
    close($stat_fh);

}

sub TAG {
    my($f, $tag) = @_;
    return "" unless $f->has_tag($tag);

    my (@values) = ($f->has_tag($tag)) ? $f->get_tag_values($tag) : ("");
    for (@values){
        #s/,/%2C/g;
    }
    my $value = join(",", @values);

#$value =~ s/^\s+(\S.*?)\s+$/$1/;

    return $value;
}



sub output_fasta{
    my ($seqHash, $cds_aa_seqs, $datadir) = @_;

    my $ffn_fh = Bio::SeqIO->new(-file=>">$datadir/cds.ffn", -format=>'fasta');
    my $faa_fh = Bio::SeqIO->new(-file=>">$datadir/cds.faa", -format=>'fasta');
    my $crispr_fh = Bio::SeqIO->new(-file=>">$datadir/crispr.ffn", -format=>'fasta');
    my $tRNA_fh = Bio::SeqIO->new(-file=>">$datadir/tRNA.ffn", -format=>'fasta');
    my $s5_fh = Bio::SeqIO->new(-file=>">$datadir/5SrRNA.ffn", -format=>'fasta');
    my $s16_fh = Bio::SeqIO->new(-file=>">$datadir/16SrRNA.ffn", -format=>'fasta');
    my $s18_fh = Bio::SeqIO->new(-file=>">$datadir/18SrRNA.ffn", -format=>'fasta');
    my $s23_fh = Bio::SeqIO->new(-file=>">$datadir/23SrRNA.ffn", -format=>'fasta');
    my $s28_fh = Bio::SeqIO->new(-file=>">$datadir/28SrRNA.ffn", -format=>'fasta');

    for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %$seqHash) {
        #for my $sid (@seqArray) {

        #my $ctg = $seqHash->{$sid}{DNA};
	next if not exists $seqHash->{$sid}{FEATURE};
        for my $f ( sort { $a->start <=> $b->start } @{ $seqHash->{$sid}{FEATURE}}) {
	    
            next if $f->primary_tag =~ /transmembrane_helix|signal_peptide/;
            my $s = $f->location->start;
	    my $e = $f->location->end;
	    my $cid = TAG($f, 'ID');
	    #msg("feature id=$cid, start=$s, end=$e");
	    my $p = $seqHash{$sid}{DNA}->trunc($f->location);
	    $p->display_id(TAG($f, 'ID') );
            my $start = $f->start;
            my $end = $f->end;
            my $strand = $f->strand;
            my $len = $f->length;
	    my $desc = "";
	    my $score = $f->score;
	    
            if ($f->primary_tag eq 'CDS'){
		if($f->has_tag('sprot_desc')){
		    $desc =  $f->has_tag('sprot_desc') ? TAG($f, 'sprot_desc') : "";
		}
		elsif($f->has_tag('tigrfam_desc')){
		    $desc =  $f->has_tag('tigrfam_desc') ? TAG($f, 'tigrfam_desc') : "";
		}
		elsif($f->has_tag('pfam_desc')){
		    $desc =  $f->has_tag('pfam_desc') ? TAG($f, 'pfam_desc') : "";
		}
		else{
		    $desc = "hypothetical protein";
		}

		my $lineage = $f->has_tag('genomedb_OC') ? TAG($f, 'genomedb_OC') : "";
		my $os = "UNKNOWN";

		if($lineage ne ""){
		    my @a = split(/;/, $lineage);
		    $os = $a[-1];
		    $os =~ s/\s+/_/g;
		}
		$desc .= " OS=$os";


		$p->desc("$desc len=$len cid=$sid") if defined $p;
                $ffn_fh->write_seq($p);

		my $aa_seq = $cds_aa_seqs->{TAG($f, 'ID')};
		$aa_seq->desc("$desc len=" . $aa_seq->length . " cid=$sid");

		$faa_fh->write_seq($aa_seq);
            }
            if ($f->primary_tag =~ /repeat_region/) {
		my $rptSeq = $f->has_tag("rpt_unit_seq") ? ($f->get_tag_values("rpt_unit_seq"))[0] : "";
		$p->desc("/numRepeat=$score /rpt_unit_seq=$rptSeq /cid=$sid") if defined $p;
                $crispr_fh->write_seq($p);
            }
            if ($f->primary_tag =~ /tRNA/) {
                $desc .= "/name=" . TAG($f, "Name");
		$p->desc("$desc /cid=$sid");
                $tRNA_fh->write_seq($p);
            }

            if ($f->primary_tag =~ /18SrRNA/) {
                $desc .= " \/Name=";
                $desc .=  $f->has_tag('Name') ? TAG($f, 'Name') : "";
                $desc .= " \/rRNA_taxon=";
                $desc .=  $f->has_tag('rRNA_taxon') ? TAG($f, 'rRNA_taxon') : "";
		
                $p->desc("$desc /cid=$sid");
                $s18_fh->write_seq($p);
            }
            if ($f->primary_tag =~ /16SrRNA/) {
                $desc .= " \/Name=";
                $desc .=  $f->has_tag('Name') ? TAG($f, 'Name') : "";
                $desc .= " \/rRNA_taxon=";
                $desc .=  $f->has_tag('rRNA_taxon') ? TAG($f, 'rRNA_taxon') : "";
		$p->desc("$desc /cid=$sid");
                $s16_fh->write_seq($p);
            }

            if ($f->primary_tag =~ /23SrRNA/) {
                $desc .= " \/Name=";
                $desc .=  $f->has_tag('Name') ? TAG($f, 'Name') : "";
                $desc .= " \/rRNA_taxon=";
                $desc .=  $f->has_tag('rRNA_taxon') ? TAG($f, 'rRNA_taxon') : "";
		$p->desc("$desc /cid=$sid");
                $p->desc($desc);
                $s23_fh->write_seq($p);
            }
	    if ($f->primary_tag =~ /28SrRNA/) {
                $desc .= " \/Name=";
                $desc .=  $f->has_tag('Name') ? TAG($f, 'Name') : "";
                $desc .= " \/rRNA_taxon=";
                $desc .=  $f->has_tag('rRNA_taxon') ? TAG($f, 'rRNA_taxon') : "";
		$p->desc("$desc /cid=$sid");
                $p->desc($desc);
                $s28_fh->write_seq($p);
            }
            if ($f->primary_tag =~ /5SrRNA|5\.8SrRNA/) {
                $desc .= " \/Name=";
                $desc .=  $f->has_tag('Name') ? TAG($f, 'Name') : "";
                $desc .= " \/rRNA_taxon=";
                $desc .=  $f->has_tag('rRNA_taxon') ? TAG($f, 'rRNA_taxon') : "";
		$p->desc("$desc /cid=$sid");
                $p->desc($desc);
                $s5_fh->write_seq($p);
            }

        }

    }

}

#use , to seperate multiple values, so, comma were escaped using %2C
sub output_tbl{
    my ($seqHash,$datadir) = @_;

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    # Write it all out!

    msg("Writing tbl output to $datadir/");
    open my $tbl_fh, '>', "$datadir/master.tbl.txt";
    for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %$seqHash){

	print $tbl_fh ">Feature $sid\n";
	next if not exists $seqHash->{$sid}{FEATURE};
	for my $f ( sort { $a->start <=> $b->start } @{ $seqHash->{$sid}{FEATURE} }) {

	    #for my $fid ( sort {$seqHash->{$sid}{$a}->start <=> $seqHash->{$sid}{$b}->start} keys %{$seqHash->{$sid}}){
	    #my $f = $seqHash->{$sid}{$fid};

	    my ($L,$R) = ($f->strand >= 0) ? ($f->start,$f->end) : ($f->end,$f->start);
	    print {$tbl_fh} "$L\t$R\t",$f->primary_tag, "\n";
	    print {$tbl_fh} "\t\t\tID\t", ($f->get_tag_values("ID"))[0], "\n" if $f->has_tag("ID");

	    foreach my  $ftag ($f->all_tags()){
		if($ftag ne "ID" 
		   && $ftag ne "score" 
		   && $ftag ne "samplePlevel" 
		   && $ftag ne "sampleDepth"
		   && $ftag ne "pfam_go" && $ftag ne "tigrfam_ec"
		   && $ftag ne "sprot_kegg" && $ftag ne "sprot_kos" && $ftag ne "sprot_ecs" && $ftag ne "sprot_go" && $ftag ne "sprot_ec"
		   && $ftag ne "foam_ecs" && $ftag ne "tigrfam_ecs" && $ftag ne "foam_acc" && $ftag ne "foam_target" && $ftag ne "foam_kos"){
		    
		    print {$tbl_fh} "\t\t\t$ftag";

		    foreach my $value ($f->get_tag_values($ftag)){
			$value ||= "";
			print {$tbl_fh} "\t$value;";
		    }
		    print {$tbl_fh} "\n";
		}
	    }

	}
    }

}
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub output_geneAnnotation{
    
    my($seqHash, $outdir) = @_;
    
    my $types = "tRNA:crispr:rRNA:sprot:tigrfam:pfam:casgene:metabolic:genomedb:ec:ko:go";
    
    my @annots = split(/:/, $types);
    my %annot_tab_file_handler = ();
    my %annot_json_file_handler = ();
    
    my @element_header = ("Gene", "Contig", "Start", "End", "Gene_Length", "Strand", "Feature");
    my @common_header = ("Gene", "Contig", "Start", "End", "Gene_Length", "Strand", "SP", "TM");
    
    my $tRNA_datatable_str = "data = [\n";
    my $rRNA_datatable_str = "data = [\n";
    my $crispr_datatable_str = "data = [\n";
    my %ec2genes = ();
    my %ko2genes = ();
    
    foreach my $type (@annots){
	print STDERR "type=$type\n";
	
	if($type =~ /tRNA|crispr|rRNA/){
	    open $annot_tab_file_handler{$type}, ">", "$outdir/$type.tab.txt" or die "Could not open $outdir/$type.tab.txt to write, $!\n";
	    open $annot_json_file_handler{$type}, ">", "$outdir/$type.datatable.json" or die "Could not open $outdir/$type.datatable.json to write, $!\n";
	}
	else{
	    open $annot_tab_file_handler{$type}, ">", "$outdir/cds.gene2$type.tab.txt" or die "Could not open $outdir/cds.gene2$type.tab.txt to write, $!\n";
	    if($type eq "ko"){
		open $annot_tab_file_handler{"gene2ko"}, ">", "$outdir/cds.gene2ko.mapping.txt" or die "Could not open $outdir/cds.gene2ko.mapping.txt to write, $!\n";
	    }
	    if($type eq "ec"){
		open $annot_tab_file_handler{"gene2ec"}, ">", "$outdir/cds.gene2ec.mapping.txt" or die "Could not open $outdir/cds.gene2ec.mapping.txt to write, $!\n";
	    }
	}
	
	
	my @header = ();
	if($type eq "crispr"){
	    push(@header, "numRepeat");
	    push(@header, "rpt_uniq_seq");
	    print {$annot_tab_file_handler{"crispr"}} "#", join("\t", (@element_header, @header)), "\n";
	}
	elsif($type eq "tRNA"){
	    push(@header, "tRNA-isoacceptor-anticodon");
	    print {$annot_tab_file_handler{"tRNA"}} "#", join("\t", (@element_header, @header)), "\n";
	    
	}
	elsif($type eq "rRNA"){
	push(@header, "Taxonomy");
	print {$annot_tab_file_handler{"rRNA"}} "#", join("\t", (@element_header, @header)), "\n";
	}
	elsif($type eq "sprot"){
	    
	    push(@header, "sprot_id");
	    push(@header, "sprot_desc");
	    push(@header, "Taxonomy");
	    print {$annot_tab_file_handler{"sprot"}} "#", join("\t", (@common_header, @header)), "\n";
	}
	elsif($type eq "pfam"){

	    push(@header, "pfam_name");
	    push(@header, "pfam_acc");
	    push(@header, "pfam_desc");
	    push(@header, "Taxonomy");
	    print {$annot_tab_file_handler{"pfam"}} "#", join("\t", (@common_header, @header)), "\n";
	}
	elsif($type eq "tigrfam"){
	    
	    push(@header, "tigrfam_acc");
	    push(@header, "tigrfam_desc");
	    push(@header, "Taxonomy");
	    print {$annot_tab_file_handler{"tigrfam"}} "#", join("\t", (@common_header, @header)), "\n";
	}
	elsif($type eq "casgene"){
	    
	    push(@header, "casgene_name");
	    push(@header, "casgene_acc");
	    push(@header, "Taxonomy");
	    print {$annot_tab_file_handler{"casgene"}} "#", join("\t", (@common_header, @header)), "\n";
	}
	elsif($type eq "metabolic"){
	    
	    push(@header, "metabolic_acc");
	    push(@header, "metabolic_process");
	    push(@header, "Taxonomy");
	    print {$annot_tab_file_handler{"metabolic"}} "#", join("\t", (@common_header, @header)), "\n";
	}
	elsif($type eq "genomedb"){
	    
	    push(@header, "genomedb_acc");
	    push(@header, "Taxonomy");
	    print {$annot_tab_file_handler{"genomedb"}} "#", join("\t", (@common_header, @header)), "\n";
	}
	elsif($type eq "ec"){
	    
	    push(@header, "ec_id");
	    push(@header, "ec_name");
	    push(@header, "Taxonomy");
	    print {$annot_tab_file_handler{"ec"}} "#", join("\t", (@common_header, @header)), "\n";
	}
	elsif($type eq "ko"){
	    
	    push(@header, "koid");
	    push(@header, "ko_name");
	    push(@header, "ko_definition");
	    push(@header, "Taxonomy");
	    print {$annot_tab_file_handler{"ko"}} "#", join("\t", (@common_header, @header)), "\n";
	    
	}
	elsif($type eq "go"){
	    
	    push(@header, "go_id");
	    push(@header, "Taxonomy");
	    print {$annot_tab_file_handler{"go"}} "#", join("\t", (@common_header, @header)), "\n";
	    
	}
    }
    
    
    for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %$seqHash) {
	next if not exists $seqHash->{$sid}{FEATURE};
        for my $f ( sort { $a->start <=> $b->start } @{ $seqHash->{$sid}{FEATURE}}) {
	    
	    my $type = $f->primary_tag;
		#next if $type !~ /$feature/;
		#next if !$f->has_tag($types);
		my $source = $f->source_tag;
		my $start = $f->start;
		my $end = $f->end;
		my $len = $end - $start + 1;
		my $strand = $f->strand;
		my $seqid = $f->seq_id;
		my $featureid = ($f->get_tag_values("ID"))[0];
		my $sp = $f->has_tag("sp") ?  ($f->get_tag_values("sp"))[0] : "No";
		my $tm = $f->has_tag("tm_num") ?  ($f->get_tag_values("tm_num"))[0] : "0";
		my $oc = $f->has_tag("genomedb_OC") ? ($f->get_tag_values("genomedb_OC"))[0] : "unknow";
		my $score = $f->score;
		my $taxon = $oc;
		$taxon =~ s/\n//;
		$taxon =~ s/\r//;
		#my $taxon = "unknow";
		#if($f->has_tag("genomedb_OC")){
		#my $value = ($f->get_tag_values("genomedb_OC"))[0];
		#($taxon) = $value =~ /(^[^;]*;[^;]*)/;
		
		#}
		my $datatable_str = "";
		#@element_header = ("Gene", "Contig", "Start", "End", "Gene_Length", "Strand", "Feature");
		
		if($type =~ /tRNA|rRNA|repeat_region/){
		    $datatable_str .= (' ' x 2) . "{\n";
		    $datatable_str .= (' ' x 2) . "\"geneid\": \"$featureid\",\n";
		    $datatable_str .= (' ' x 2) . "\"contigid\": \"$seqid\",\n";
		    $datatable_str .= (' ' x 2) . "\"start\": \"$start\",\n";
		    $datatable_str .= (' ' x 2) . "\"end\": \"$end\",\n";
		    $datatable_str .= (' ' x 2) . "\"genelen\": \"$len\",\n";
		    $datatable_str .= (' ' x 2) . "\"strand\": \"$strand\",\n";
		    $datatable_str .= (' ' x 2) . "\"feature\": \"$type\",\n";
		}
		
		
		if(exists $annot_tab_file_handler{"crispr"} && $type eq "repeat_region"){
		    my $count = $f->get_tag_values("rpt_unit_seq");
		    for (my $i = 0; $i < $count; $i++){
			my $rpt_unit_seq = ($f->get_tag_values("rpt_unit_seq"))[$i];
			my @row = ($featureid, $seqid, $start, $end, $len, $strand,$type, $score,$rpt_unit_seq);
			print {$annot_tab_file_handler{"crispr"}} join("\t", @row), "\n";
			$crispr_datatable_str .= $datatable_str;
			$crispr_datatable_str .= (' ' x 2) . "\"numRepeat\": \"$score\",\n";
			$crispr_datatable_str .= (' ' x 2) . "\"rpt_unit_seq\": \"$rpt_unit_seq\"\n";
			$crispr_datatable_str .=(' ' x 2) . "},\n";
		    }
		    
		}
		if(exists $annot_tab_file_handler{"tRNA"} && $type eq "tRNA"){
		    my $count = $f->has_tag("name") ? $f->get_tag_values("name") : 0;
		    for (my $i = 0; $i < $count; $i++){
			my $name = ($f->get_tag_values("name"))[$i];
			my @row = ($featureid, $seqid, $start, $end, $len, $strand,$type, $name);
			print {$annot_tab_file_handler{"tRNA"}} join("\t", @row), "\n";
			
			$tRNA_datatable_str .= $datatable_str;
			$tRNA_datatable_str .= (' ' x 2) . "\"name\": \"$name\"\n";
			$tRNA_datatable_str .=(' ' x 2) . "},\n";
		    }
		    
		}
		if(exists $annot_tab_file_handler{"rRNA"} && $type =~ /SrRNA/){
		    my $count = $f->get_tag_values("rRNA_taxon");
		    for (my $i = 0; $i < $count; $i++){
			#my $name = ($f->get_tag_values("Name"))[$i];
			my $taxon = ($f->get_tag_values("rRNA_taxon"))[$i];
			my @row = ($featureid, $seqid, $start, $end, $len, $strand,$type, $taxon);
			print {$annot_tab_file_handler{"rRNA"}} join("\t", @row), "\n";
			$rRNA_datatable_str .= $datatable_str;
			$rRNA_datatable_str .= (' ' x 2) . "\"taxon\": \"$taxon\"\n";
			$rRNA_datatable_str .=(' ' x 2) . "},\n";
			
		    }
		    
		}
		if(exists $annot_tab_file_handler{"sprot"} && $f->has_tag("sprot_target")){
		    my $count = $f->get_tag_values("sprot_id");
		    for (my $i = 0; $i < $count; $i++){
			my $sprotid = ($f->get_tag_values("sprot_id"))[$i];
			my $desc = ($f->get_tag_values("sprot_desc"))[$i];
			$desc =~ s/\n//;
			$taxon =~ s/\n//g;
			my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $sprotid, "\"$desc\"", $taxon);
			print {$annot_tab_file_handler{"sprot"}} join("\t", @row), "\n";
			
		    }
		    
		}
		if(exists $annot_tab_file_handler{"pfam"} && $f->has_tag("pfam_acc")){
		    my $count = $f->get_tag_values("pfam_acc");
		    for (my $i = 0; $i < $count; $i++){
			my $acc = ($f->get_tag_values("pfam_acc"))[$i];
			my $pid = ($f->get_tag_values("pfam_id"))[$i];
			my $desc = "\"" . ($f->get_tag_values("pfam_desc"))[$i] . "\"";
			
			my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $pid, $acc, $desc,$taxon);
			print {$annot_tab_file_handler{"pfam"}} join("\t", @row), "\n";
			
		    }
		}
		if(exists $annot_tab_file_handler{"tigrfam"}  && $f->has_tag("tigrfam_acc")){
		    my $count = $f->get_tag_values("tigrfam_acc");
		    for (my $i = 0; $i < $count; $i++){
			my $acc = ($f->get_tag_values("tigrfam_acc"))[$i];
			my $desc = "\"" . ($f->get_tag_values("tigrfam_desc"))[$i] . "\"";
			my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $acc, $desc,$taxon);
			print {$annot_tab_file_handler{"tigrfam"}} join("\t", @row), "\n";
		    }
		    
		}
		if(exists $annot_tab_file_handler{"casgene"} && $f->has_tag("casgene_acc")){
		    
		    my $count = $f->get_tag_values("casgene_acc");
		    for (my $i = 0; $i < $count; $i++){
			my $name = ($f->get_tag_values("casgene_name"))[$i];
			my $acc = ($f->get_tag_values("casgene_acc"))[$i];
			$acc =~ s/\r//;
			$acc =~ s/\n//;
			my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $name, "\"$acc\"",$taxon);
			print {$annot_tab_file_handler{"casgene"}} join("\t", @row), "\n";
		    }
		    
		    
		}
		if(exists $annot_tab_file_handler{"metabolic"} && $f->has_tag("metabolic_acc")){
		    
		    my $count = $f->get_tag_values("metabolic_acc");
		    for (my $i = 0; $i < $count; $i++){
			my $acc = ($f->get_tag_values("metabolic_acc"))[$i];
			my $process = "\"" . ($f->get_tag_values("metabolic_process"))[$i] . "\"";
			my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $acc, $process,$taxon);
			print {$annot_tab_file_handler{"metabolic"}} join("\t", @row), "\n";
		    }
		    
		    
		}
		if(exists $annot_tab_file_handler{"genomedb"} && $f->has_tag("genomedb_acc")){
		    my $count = $f->get_tag_values("genomedb_acc");
		    for (my $i = 0; $i < $count; $i++){
			my $acc = ($f->get_tag_values("genomedb_acc"))[$i];
			my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $acc,$taxon);
			print {$annot_tab_file_handler{"genomedb"}} join("\t", @row), "\n";
		    }
		}
		
		if(exists $annot_tab_file_handler{"ec"} && $f->has_tag("allec_ids")){
		    my $count = $f->get_tag_values("allec_ids");
		    for (my $i = 0; $i < $count; $i++){
			my $ec = ($f->get_tag_values("allec_ids"))[$i];
			my $de = "";
			my $stmt = "select de from enzyme where ec=\"$ec\"";
			my $sth = $dbh->prepare( $stmt);
			$sth->execute();
			my $all = $sth->fetchall_arrayref();
			foreach my $row (@$all) {
			    ($de) = @$row;
			}
			$sth->finish();
			my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $ec,"\"$de\"","\"$taxon\"");
			print {$annot_tab_file_handler{"ec"}} join("\t", @row), "\n";
			print {$annot_tab_file_handler{"gene2ec"}}  "$featureid\t$ec\n";
			$ec2genes{$ec}->{$featureid} = 1;
		    }
		}
		
		if(exists $annot_tab_file_handler{"ko"} && $f->has_tag("allko_ids")){
		    my $count = $f->get_tag_values("allko_ids");
		    for (my $i = 0; $i < $count; $i++){
			my $ko = ($f->get_tag_values("allko_ids"))[$i];
			my $koname = "";
			my $kode = "";
			
			my $stmt = "select name, de from kos where koid=\"$ko\"";
			my $sth = $dbh->prepare( $stmt);
			$sth->execute();
			my $all = $sth->fetchall_arrayref();
			
			foreach my $row (@$all) {
			    ($koname, $kode) = @$row;
			}
			$sth->finish();
			


			my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $ko,"\"$koname\"", "\"$kode\"","\"$taxon\"");
			print {$annot_tab_file_handler{"ko"}} join("\t", @row), "\n";
			print {$annot_tab_file_handler{"gene2ko"}}  "$featureid\t$ko\n";
			$ko2genes{$ko}->{$featureid} = 1;
			
		    }
		}
		if(exists $annot_tab_file_handler{"go"} && $f->has_tag("allgo_ids")){
		    my $count = $f->get_tag_values("allgo_ids");
		    for (my $i = 0; $i < $count; $i++){
			my $go = ($f->get_tag_values("allgo_ids"))[$i];
			my $name = "";
			my $stmt = "select name from Gos where goid=\"$go\"";
			my $sth = $dbh->prepare( $stmt);
			$sth->execute();
			my $all = $sth->fetchall_arrayref();
			
			foreach my $row (@$all) {
			    ($name) = @$row;
			}
			$sth->finish();
			

			my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $go, "\"$name\"", "\"$taxon\"");
			print {$annot_tab_file_handler{"go"}} join("\t", @row), "\n";
		
		    }
		}
	    }
	
    }
    $crispr_datatable_str =~ s/,$//;
    $crispr_datatable_str .= (' ' x 1) . "]\n";
    print {$annot_json_file_handler{"crispr"}} $crispr_datatable_str;
    close($annot_json_file_handler{"crispr"});
    $tRNA_datatable_str =~ s/,$//;
    $tRNA_datatable_str .= (' ' x 1) . "]\n";
    print {$annot_json_file_handler{"tRNA"}} $tRNA_datatable_str;
    close($annot_json_file_handler{"tRNA"});
    $rRNA_datatable_str =~ s/,$//;
    $rRNA_datatable_str .= (' ' x 1) . "]\n";
    print {$annot_json_file_handler{"rRNA"}} $rRNA_datatable_str;
    close($annot_json_file_handler{"rRNA"});
    
    map{close($annot_tab_file_handler{$_})} keys %annot_tab_file_handler;
    return (\%ko2genes, \%ec2genes);
}
sub predict_kegg_pathways{
    
    my ($prefix, $ko2genes, $seqHash) = @_;
    my $cmd = "MinPath1.4.py -ko $prefix.mapping.txt -report $prefix.minpath -details $prefix.minpath.details > /dev/null 2>&1;";
    msg("******start running minpath $cmd\n");
    runcmd("$cmd");
    
    open(INPUT, "$prefix.minpath.details") or die "Could not find $prefix.minpath.details to read\n";
    $/= "\npath";
    
   #path 00010 fam0 56 fam-found 36 # Glycolysis / Gluconeogenesis
    #K00001 hits 22 # E1.1.1.1, adh
    
    #my %gene2pathways = ();
    my %pathways = ();
    
    while (<INPUT>) {
	my $fams = "";
	my @items = split(/\n/, $_);
	my $geneCount = 0;
	my $pid = "";
	foreach my $item (@items){
	    $item =~ s/^\s+//;
	    if(my($id, $total_family_num, $total_family_found, $name) = $item =~ /(\d+)\s+fam0\s+(\d+)\s+fam-found\s+(\d+)\s+\#\s+(\S.*?)$/){
		$pid = $id;
		$pathways{$pid}->{fam_all} = $total_family_num;
		$pathways{$pid}->{fam_found} = $total_family_found;
		$pathways{$pid}->{pname} = $name;
		
	    }
	    elsif(my($ko, $geneCount, $koname) = $item =~ /^(\S+)\s+hits\s+(\d+?)\s+\#?\s+(\S.*)$/){
		for my $geneid (keys %{$ko2genes->{$ko}}){
		    $gene2pathways{$geneid}->{KEGG}->{$pid}++;
		    $pathways{$pid}->{$ko}->{hits}=$geneCount;
		    $pathways{$pid}->{$ko}->{name}=$koname;
		}
	    }
	}
    }
    close(INPUT);
    $/= "\n";
    open(GENE2PATHWAY, ">$datadir/cds.gene2kegg.tab.txt") or die "Could not open $datadir/cds.gene2kegg.tab.txt to write, $!\n";
    print GENE2PATHWAY "#Geneid\tKegg_id\tPathway_name\tTaxonomy\n";
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    my $geneid = ($f->get_tag_values("ID"))[0];
	    
	    next unless $f->primary_tag eq "CDS";
	    my $oc = $f->has_tag("genomedb_OC") ? ($f->get_tag_values("genomedb_OC"))[0] : "unknow";
	    if(exists $gene2pathways{$geneid}->{KEGG}){
		
		foreach my $pid (keys %{$gene2pathways{$geneid}->{KEGG}}){
		    $f->add_tag_value("kegg_pathway_id", $pid);
		    $f->add_tag_value("kegg_pathway_name", $pathways{$pid}->{pname});
		    
		    print GENE2PATHWAY "$geneid\t$pid\t\"$pathways{$pid}->{pname}\"\t\"$oc\"\n";
		    if($gene2pathways{$geneid}->{KEGG}->{$pid} > 1){

			#print "********$geneid $pid $gene2pathways{$geneid}->{$pid} times\n";
		    }
		}
		
	    }
	}
    }
    close(GENE2PATHWAY);
    
    foreach my $pid (keys %pathways){
	my $tn = exists $pathways{$pid}->{fam_all} ? $pathways{$pid}->{fam_all} : 0;
	my $tfn = exists $pathways{$pid}->{fam_found} ? $pathways{$pid}->{fam_found} : 0;
	my $pname = exists $pathways{$pid}->{pname} ? $pathways{$pid}->{pname} : "";
	my @kos = grep !/fam_all|fam_found|pname/, sort keys %{$pathways{$pid}};
	
	foreach my $ko (@kos){
	    my $hits = $pathways{$pid}->{$ko}->{hits};
	    my $koname = $pathways{$pid}->{$ko}->{name};
	 
	}
	
    }
    
    return \%pathways;
}

sub predict_metacyc_pathways{
    
    my ($prefix, $ec2genes, $seqHash) = @_;

    open (DATA, "$txtdir/MetaCyc_pathways_id2name.tab.txt") or die "Could not open $txtdir/MetaCyc_pathways_id2name.tab.txt file to read\n";
    
    #UNIQUE-ID	TYPES	COMMON-NAME
    my %pid2info = ();
    while (<DATA>) {
	chomp;
		
	next if /^#/;
	next if !/\S/;
	my @line = split(/\t/, $_);
	my $id = $line[0];
	my $type = $line[1];
	my $name = $line[2];
	#$type =~ s/<i>|<\/i>//g;
	#$name =~ s/<i>|<\/i>//g;
	$pid2info{$id}->{ptype} = $line[1];
	$pid2info{$id}->{pname} = $line[2];
	
    }
    close(DATA);
    


    my $cmd = "MinPath1.4.py -any $prefix.mapping.txt -map ec2path -report $prefix.minpath -details $prefix.minpath.details > /dev/null 2>&1;";
    msg("******start running minpath $cmd\n");
    runcmd("$cmd");
    
    open(INPUT, "$prefix.minpath.details") or die "Could not find $prefix.minpath.details to read\n";
    $/= "\npath";
    
   #path 00010 fam0 56 fam-found 36 # Glycolysis / Gluconeogenesis
    #K00001 hits 22 # E1.1.1.1, adh
    
    #my %gene2pathways = ();
    my %pathways = ();
    
    while (<INPUT>) {
	my $fams = "";
	my @items = split(/\n/, $_);
	my $geneCount = 0;
	my $pid = "";
	foreach my $item (@items){
	    $item =~ s/^\s+//;
	    if(my($total_family_num, $total_family_found, $id) = $item =~ /fam0\s+(\d+)\s+fam-found\s+(\d+)\s+\#\s+(\S.*?)$/){
		$pid = $id;
		$pathways{$pid}->{fam_all} = $total_family_num;
		$pathways{$pid}->{fam_found} = $total_family_found;
		
	    }
	    elsif(my($ec, $geneCount) = $item =~ /^(\S+)\s+hits\s+(\d+?)/){
		for my $geneid (keys %{$ec2genes->{$ec}}){
		    $gene2pathways{$geneid}->{metacyc}->{$pid}++;
		    $pathways{$pid}->{$ec}->{hits}=$geneCount;
		}
	    }
	}
    }
    close(INPUT);
    $/= "\npath";
    open(GENE2PATHWAY, ">$datadir/cds.gene2metacyc.tab.txt") or die "Could not open $datadir/cds.gene2metacyc.tab.txt to write, $!\n";
    print GENE2PATHWAY "#Geneid\tMetaCyc_id\tPathway_type\tPathway_name\n";
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    my $geneid = ($f->get_tag_values("ID"))[0];
	    
	    next unless $f->primary_tag eq "CDS";
	    my $oc = $f->has_tag("genomedb_OC") ? ($f->get_tag_values("genomedb_OC"))[0] : "unknow";
	    if(exists $gene2pathways{$geneid}->{metacyc}){
		
		foreach my $pid (keys %{$gene2pathways{$geneid}->{metacyc}}){
		    #print "*****$pid***\t$pid2info{$pid}->{types}\t$pid2info{$pid}->{name}\n";
		    $f->add_tag_value("metacyc_pathway_id", $pid);
		    my $ptype = exists $pid2info{$pid}->{ptype} ? $pid2info{$pid}->{ptype} : "";
		    $f->add_tag_value("metacyc_pathway_type", $ptype);
		    my $pname = exists $pid2info{$pid}->{pname} ? $pid2info{$pid}->{pname} : "";
		    $f->add_tag_value("metacyc_pathway_name", $pname);
		    print GENE2PATHWAY "$geneid\t$pid\t\"$ptype\"\t\"$pname\"\t\"$oc\"\n";
		}
		
	    }
	}
    }
    close(GENE2PATHWAY);
  
    foreach my $pid (keys %pathways){
	my $tn = $pathways{$pid}->{fam_all};
	my $tfn = $pathways{$pid}->{fam_found};
	my $pname = exists $pid2info{$pid}->{pname} ? $pid2info{$pid}->{pname} : "";
	my $ptype = exists $pid2info{$pid}->{ptype} ? $pid2info{$pid}->{ptype} : "";
	$pathways{$pid}->{ptype} = $ptype;
	$pathways{$pid}->{pname} = $pname;
	
	my @ecs = grep !/fam_all|fam_found|ptype|pname/, sort keys %{$pathways{$pid}};
	foreach my $ec (@ecs){
	    my $hits = $pathways{$pid}->{$ec}->{hits};
	    
	}
	
    }
    
    return \%pathways;
    
}

sub output_profiles{

    my($seqHash, $outdir) = @_;
    
    my $bin = "$FindBin::RealBin";
    my %taxon2genes_cds = ();
    my %taxon2genes_ssu = ();
    my %taxon2genes_lsu = ();
    my %pfam2genes = ();
    
    my %pfam = ();
    my %tigrfam2genes = ();
    my %tigrfam = ();
    my %ko = ();
    my %ko2genes = ();
    my %ec = ();
    my %ec2genes = ();

    my %go = ();
    my %go2genes = ();

    my %stats = ();
    my %metabolic2genes = ();
    my %metabolic = ();
    
    my %keggs = ();
    my %keggPathway2genes = ();
    my %metacycPathway2genes = ();
    
    
    for my $sid (keys %$seqHash) {
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
    
	    my $feature = $f->primary_tag;
	    my $source = $f->source_tag;
	    my $start = $f->start;
	    my $end = $f->end;
	    my $len = $end - $start + 1;
	    my $strand = $f->strand;
	    my $sid = $f->seq_id;
	    my $geneid = ($f->get_tag_values("ID"))[0];
	    my $oc = $f->has_tag("genomedb_OC") ? ($f->get_tag_values("genomedb_OC"))[0] : "unknow";
	    my $taxon = $oc;
	    
	    my $mdepth_col_count = $f->has_tag("mdepth_cols") ? $f->get_tag_values("mdepth_cols") : -1;
	    my %sample2depth = ();
	    
	    for (my $i = 0; $i < $mdepth_col_count; $i++){
		my $name = ($f->get_tag_values("mdepth_cols"))[$i];
		my $depth = ($f->get_tag_values("mdepth_values"))[$i];
		$stats{cds}->{$name} += $depth;
		$sample2depth{$name} = $depth;
		
	    }
	
	    if($feature =~ /CDS/){
		
		$stats{cds}->{totalCount}++;
		$taxon2genes_cds{$taxon}->{count}++;
		
		foreach my $name (keys %sample2depth){
		    my $depth = $sample2depth{$name};
		    $taxon2genes_cds{$taxon}->{$name} += $depth;
		}
		
		if($f->has_tag("pfam_acc")){
		    my $acc_count = $f->get_tag_values("pfam_acc");
		    for (my $i = 0; $i < $acc_count; $i++){
			my $acc = ($f->get_tag_values("pfam_acc"))[$i];
			my $id = ($f->get_tag_values("pfam_id"))[$i];
			my $desc = ($f->get_tag_values("pfam_desc"))[$i];
			
			if(not exists $pfam{$acc}){
			    $pfam{$acc}->{ID} = $id;
			    $pfam{$acc}->{Desc} = $desc;
			}
			$pfam2genes{$acc}->{count}++;
			foreach my $name (keys %sample2depth){
			    my $depth = $sample2depth{$name};
			    $pfam2genes{$acc}->{$name} += $depth;
			}
		    }
		}
		if($f->has_tag("tigrfam_acc")){
	    
		    my $acc_count = $f->get_tag_values("tigrfam_acc");
		    for (my $i = 0; $i < $acc_count; $i++){
			my $acc = ($f->get_tag_values("tigrfam_acc"))[$i];
			my $name = ($f->get_tag_values("tigrfam_name"))[$i];
			my $desc = ($f->get_tag_values("tigrfam_desc"))[$i];
			#my $mainrole = $f->has_tag("tigrfam_mainrole") ? ($f->get_tag_values("tigrfam_mainrole"))[$i] : "";
			if(not exists $tigrfam{$acc}){
			    $tigrfam{$acc}->{Name} = $name;
			    $tigrfam{$acc}->{Desc} = $desc;
			    #$tigrfam{$acc}->{Mainrole} = $mainrole;
			}
			$tigrfam2genes{$acc}->{count}++;
			foreach my $name (keys %sample2depth){
			    my $depth = $sample2depth{$name};
			    $tigrfam2genes{$acc}->{$name} += $depth;
			}
		    }
		}
		if($f->has_tag("allko_ids")){
		    my $ko_count = $f->get_tag_values("allko_ids");
		    for (my $i = 0; $i < $ko_count; $i++){
			my $koid = ($f->get_tag_values("allko_ids"))[$i];
			my $name = "";
			my $de = "";
			my $stmt = "select name,de from kos where koid=\"$koid\"";
			my $sth = $dbh->prepare( $stmt);
			$sth->execute();
			my $all = $sth->fetchall_arrayref();
			foreach my $row (@$all) {
			    ($name, $de) = @$row;
			}
			$sth->finish();
			if(not exists $ko{$koid}){
			    $ko{$koid}->{name} = $name;
			    $ko{$koid}->{de} = $de;
			}
			$ko2genes{$koid}->{count}++;
			foreach my $name (keys %sample2depth){
			    my $depth = $sample2depth{$name};
			    $ko2genes{$koid}->{$name} += $depth;
			}
		    }
		}
		if($f->has_tag("allec_ids")){
		    my $ec_count = $f->get_tag_values("allec_ids");
		    for (my $i = 0; $i < $ec_count; $i++){
			my $ecid = ($f->get_tag_values("allec_ids"))[$i];
			my $de = "";
			my $stmt = "select de from enzyme where ec=\"$ecid\"";
			my $sth = $dbh->prepare( $stmt);
			$sth->execute();
			my $all = $sth->fetchall_arrayref();
			foreach my $row (@$all) {
			    ($de) = @$row;
			}
			$sth->finish();
			if(not exists $ec{$ecid}){
			    $ec{$ecid}->{de} = $de;
			}
			$ec2genes{$ecid}->{count}++;
			foreach my $name (keys %sample2depth){
			    my $depth = $sample2depth{$name};
			    $ec2genes{$ecid}->{$name} += $depth;
			}
		    }
		}
		if($f->has_tag("allgo_ids")){
		    my $go_count = $f->get_tag_values("allgo_ids");
		    for (my $i = 0; $i < $go_count; $i++){
			my $goid = ($f->get_tag_values("allgo_ids"))[$i];
			my $name = "";
			my $stmt = "select name from Gos where goid=\"$goid\"";
			my $sth = $dbh->prepare( $stmt);
			$sth->execute();
			my $all = $sth->fetchall_arrayref();
			foreach my $row (@$all) {
			    ($name) = @$row;
			}
			$sth->finish();
			
			if(not exists $go{$goid}){
			    $go{$goid}->{name} = $name;
			}
			$go2genes{$goid}->{count}++;
			foreach my $name (keys %sample2depth){
			    my $depth = $sample2depth{$name};
			    $go2genes{$goid}->{$name} += $depth;
			}
		    }
		}
		if($f->has_tag("metabolic_acc")){
		    
		    my $acc_count = $f->get_tag_values("metabolic_acc");
		    for (my $i = 0; $i < $acc_count; $i++){
			my $acc = ($f->get_tag_values("metabolic_acc"))[$i];
			#compound:Sulfur;process:thiosulfate oxidation;gene:soxY;
			my $desc = ($f->get_tag_values("metabolic_process"))[$i];
			my($compound, $proc) = $desc =~ /^compound:(\S.*?);process:(\S+.*?).gene:/;
			$desc =~ s/compound:|process:|gene://g;
			if(not exists $metabolic{$acc}){
			    $metabolic{$acc}->{Process} = $proc;
			    $metabolic{$acc}->{Compound} = $compound;
			    $metabolic{$acc}->{Desc} = $desc;
			}
			
			$metabolic2genes{$acc}->{count}++;
			foreach my $name (keys %sample2depth){
			    my $depth = $sample2depth{$name};
			    $metabolic2genes{$acc}->{$name} += $depth;
			}
		    }
		}
		if($f->has_tag("kegg_pathway_id")){
		    my $pathway_count = $f->get_tag_values("kegg_pathway_id");
		    for (my $i = 0; $i < $pathway_count; $i++){
			my $pid = ($f->get_tag_values("kegg_pathway_id"))[$i];
			my $pname = ($f->get_tag_values("kegg_pathway_name"))[$i];
			$keggPathway2genes{$pid}->{count}++;
			#$keggPathway2genes{$pid}->{pname} = $pname;
			foreach my $name (keys %sample2depth){
			    my $depth = $sample2depth{$name};
			    $keggPathway2genes{$pid}->{$name} += $depth;
			}
		    }
		}
		if($f->has_tag("metacyc_pathway_id")){
		    my $pathway_count = $f->get_tag_values("metacyc_pathway_id");
		    for (my $i = 0; $i < $pathway_count; $i++){
			my $pid = ($f->get_tag_values("metacyc_pathway_id"))[$i];
			$metacycPathway2genes{$pid}->{count}++;
			foreach my $name (keys %sample2depth){
			    my $depth = $sample2depth{$name};
			    $metacycPathway2genes{$pid}->{$name} += $depth;
			}
		    }
		}
	    }
	    
	    elsif($feature =~ /16SrRNA|18SrRNA/){
		my $taxon = ($f->get_tag_values("rRNA_taxon"))[0];
		$stats{ssu}->{totalCount}++;
		$taxon2genes_ssu{$taxon}->{count}++;
		
		foreach my $name (keys %sample2depth){
		    my $depth = $sample2depth{$name};
		    $taxon2genes_ssu{$taxon}->{$name} += $depth;
		    $stats{ssu}->{$name} += $depth;
		}
	    }
	    elsif($feature =~ /23SrRNA|28SrRNA/){
		
		my $taxon = ($f->get_tag_values("rRNA_taxon"))[0];
		$stats{lsu}->{totalCount}++;
		$taxon2genes_lsu{$taxon}->{count}++;
		
		foreach my $name (keys %sample2depth){
		    my $depth = $sample2depth{$name};
		    $taxon2genes_lsu{$taxon}->{$name} += $depth;
		    $stats{lsu}->{$name} += $depth;
		    
		}
	    }
	}
    }
    
    output_profile(\%stats, \%taxon2genes_cds, "cds", "taxon", $outdir);
    output_profile(\%stats, \%taxon2genes_ssu, "ssu", "taxon", $outdir, );
    output_profile(\%stats, \%taxon2genes_lsu, "lsu", "taxon", $outdir);
    output_profile(\%stats, \%pfam2genes, "cds", "pfam", $outdir, \%pfam);
    output_profile(\%stats, \%tigrfam2genes, "cds", "tigrfam", $outdir, \%tigrfam);
    output_profile(\%stats, \%metabolic2genes, "cds", "metabolic", $outdir, \%metabolic);
    output_profile(\%stats, \%ko2genes, "cds", "ko", $outdir, \%ko);
    output_profile(\%stats, \%ec2genes, "cds", "ec", $outdir, \%ec);
    output_profile(\%stats, \%go2genes, "cds", "go", $outdir, \%go);

    output_pathway_profile(\%stats, \%keggPathway2genes, "cds", "kegg", $outdir, $keggs);
    output_pathway_profile(\%stats, \%metacycPathway2genes, "cds", "metacyc", $outdir, $metacycs);
    
    my @msamples = grep !/totalCount|totalAvgDepth/, sort keys %{$stats{ssu}};
    
###########SSU
    my $cmd = "$^X $bin/output_tree_json.pl -i $outdir/taxon.ssu.profile.tab.txt -t Taxonomy> $outdir/taxon.ssu.tree.json;";
    $cmd .= "$^X $bin/output_sunburst_json.pl -i $outdir/taxon.ssu.profile.tab.txt -t Taxonomy -p $outdir/taxon.ssu.profile.sunburst;";
    if(@msamples > 0){
	$cmd .= "$^X $bin/output_sunburst_json.abund.pl -i $outdir/taxon.ssu.profile.tab.txt -t Taxonomy -p $outdir/taxon.ssu.profile.sunburst";
    }
    msg("******start running $cmd\n");
    runcmd("$cmd");
    msg("******Finish running $cmd \n\n");
    output_treecol("$outdir/taxon.ssu.treecol.json", "Taxonomy", \@msamples);
    
##########LSU
    $cmd = "$^X $bin/output_tree_json.pl -i $outdir/taxon.lsu.profile.tab.txt -t Taxonomy> $outdir/taxon.lsu.tree.json;";
    $cmd .= "$^X $bin/output_sunburst_json.pl -i $outdir/taxon.lsu.profile.tab.txt -t Taxonomy -p $outdir/taxon.lsu.profile.sunburst;";
    if(@msamples > 0){
	$cmd .= "$^X $bin/output_sunburst_json.abund.pl -i $outdir/taxon.lsu.profile.tab.txt -t Taxonomy -p $outdir/taxon.lsu.profile.sunburst";
    }
    msg("******start running $cmd\n");
    runcmd("$cmd");
    msg("******Finish running $cmd \n\n");
    output_treecol("$outdir/taxon.lsu.treecol.json", "Taxonomy", \@msamples);
    
##########CDS
    $cmd = "$^X $bin/output_tree_json.pl -i $outdir/taxon.cds.profile.tab.txt -t Taxonomy> $outdir/taxon.cds.tree.json;";
    $cmd .= "$^X $bin/output_sunburst_json.pl -i $outdir/taxon.cds.profile.tab.txt -t Taxonomy -p $outdir/taxon.cds.profile.sunburst;";
    if(@msamples > 0){
	$cmd .= "$^X $bin/output_sunburst_json.abund.pl -i $outdir/taxon.cds.profile.tab.txt -t Taxonomy -p $outdir/taxon.cds.profile.sunburst";
    }
    
    msg("******start running $cmd\n");
    runcmd("$cmd");
    msg("******Finish running $cmd \n\n");
    output_treecol("$outdir/taxon.cds.treecol.json","Taxonomy", \@msamples);
    
    
###########metabolic                                                                                                                                                                                               
    $cmd = "$^X $bin/output_tree_json.pl -i $outdir/metabolic.cds.profile.tab.txt -t Compound-process-gene> $outdir/metabolic.cds.profile.tree.json;";
    $cmd .= "$^X $bin/output_sunburst_json.pl -i $outdir/metabolic.cds.profile.tab.txt -t Compound-process-gene -p $outdir/metabolic.cds.profile.sunburst;";
    if(@msamples > 0){
	$cmd .= "$^X $bin/output_sunburst_json.abund.pl -i $outdir/metabolic.cds.profile.tab.txt -t  Foam -p $outdir/metabolic.cds.profile.sunburst";
    }
    
    msg("******start running $cmd\n");
    runcmd("$cmd");
    msg("******Finish running $cmd \n\n");
    output_treecol("$outdir/metabolic.cds.treecol.json", "Metabolic_genes", \@msamples);
 
}

sub output_profile{
    my($sref, $dataref, $feature, $attr, $outdir, $iddesc) = @_;
    my %stats = %$sref;
    my %data = %$dataref;
    my $total_count = $stats{$feature}->{totalCount};
    my $total_depth = exists $stats{$feature}->{totalAvgDepth} ?  $stats{$feature}->{totalAvgDepth} : 0;
    my @msamples = grep !/totalCount|totalAvgDepth/, sort keys %{$stats{$feature}};
    
    my @header = ("Count", "Count_pct", "Abund", "Abund_pct");
    unshift(@header, "#Taxon") if($attr =~ /taxon/);
    unshift(@header, "#Accession", "ID", "Description") if $attr =~ /pfam/;
    unshift(@header, "#Accession", "Name", "Function") if $attr =~ /tigrfam/;
    unshift(@header, "#KO", "Name", "Definition") if $attr =~ /ko/;
    unshift(@header, "#Enzyme_id", "Name") if $attr =~ /ec/;
    unshift(@header, "#GO_id", "Name") if $attr =~ /go/;
    unshift(@header, "#Gene-Compound-Process") if $attr =~ /metabolic/;
    foreach (@msamples){
	push(@header, "Abund\_$_");
	push(@header, "Abund\_pct_$_");
    }
    
    open  my $data_tab, ">", "$outdir/$attr.$feature.profile.tab.txt" or die "Could not open $outdir/$attr.$feature.profile.tab.txt to write, $!\n";
    print $data_tab join("\t", @header), "\n";
    
    open my $data_json, ">", "$outdir/$attr.$feature.profile.datatable.json" or die "Could not open $outdir/$attr.$feature.profile.datatable.json to write, $!\n";

    #output tablecol
    output_tablecols($data_json, \@msamples, $attr);
    
    print $data_json "\ndata = [\n";
    my $keyCount = keys %data;
    
    my $index = 0;
    #my @cols = ();
    for my $key (keys %data){
	my @cols = ();
	push(@cols, $key) if $attr !~ /metabolic/;
	print $data_json (' ' x 2) . "{\n";
	
	print $data_json (' ' x 3) . "\"taxon\": \"$key\",\n" if($attr =~ /taxon/);
	print $data_json (' ' x 3) . "\"accession\": \"$key\",\n" if($attr =~ /pfam|tigrfam/);
	print $data_json (' ' x 3) . "\"gene\": \"$key\",\n" if($attr =~ /metabolic/);
	print $data_json (' ' x 3) . "\"koid\": \"$key\",\n" if($attr =~ /ko/);
	print $data_json (' ' x 3) . "\"ecid\": \"$key\",\n" if($attr =~ /ec/);
	print $data_json (' ' x 3) . "\"goid\": \"$key\",\n" if($attr =~ /go/);
	if(defined $iddesc){
	    if($attr eq "pfam"){
		my $id = exists $iddesc->{$key}->{ID} ? $iddesc->{$key}->{ID} : "";
		my $desc = exists $iddesc->{$key}->{Desc} ? $iddesc->{$key}->{Desc} : "";
		push(@cols, ($id, "\"$desc\""));
		
		print $data_json (' ' x 3) . "\"id\": \"$id\",\n";
		print $data_json (' ' x 3) . "\"desc\": \"$desc\",\n";
	    }
	    elsif($attr eq "tigrfam"){
		my $name = exists $iddesc->{$key}->{Name} ? $iddesc->{$key}->{Name} : "";
		my $desc = exists $iddesc->{$key}->{Desc} ? $iddesc->{$key}->{Desc} : "";
		push(@cols, ($name, "\"$desc\""));
		print $data_json (' ' x 3) . "\"name\": \"$name\",\n";
		print $data_json (' ' x 3) . "\"desc\": \"$desc\",\n";
	    }
	    elsif($attr eq "ko"){
		my $name = exists $iddesc->{$key}->{name} ? $iddesc->{$key}->{name} : "";
		my $de = exists $iddesc->{$key}->{de} ? $iddesc->{$key}->{de} : "";
		$name =~ s/\""/\\\"/g;
		$de =~ s/\""/\\\"/g;
		
		push(@cols, ($name, "\"$de\""));
		print $data_json (' ' x 3) . "\"name\": \"$name\",\n";
		print $data_json (' ' x 3) . "\"definition\": \"$de\",\n";
	    }
	    elsif($attr eq "ec"){
		my $de = exists $iddesc->{$key}->{de} ? $iddesc->{$key}->{de} : "";
		$de =~ s/\""/\\\"/g;
		push(@cols, ("\"$de\""));
		print $data_json (' ' x 3) . "\"definition\": \"$de\",\n";
	    }
	    elsif($attr eq "go"){
		my $name = exists $iddesc->{$key}->{name} ? $iddesc->{$key}->{name} : "";
		$name =~ s/\""/\\\"/g;
		push(@cols, ("\"$name\""));
		print $data_json (' ' x 3) . "\"name\": \"$name\",\n";
	    }
	    
	    elsif($attr eq "metabolic"){
                my $compound = exists $iddesc->{$key}->{Compound} ? $iddesc->{$key}->{Compound} : "";
                my $process = exists $iddesc->{$key}->{Process} ? $iddesc->{$key}->{Process} : "";
		my $desc = exists $iddesc->{$key}->{Desc} ? $iddesc->{$key}->{Desc} : "";
                push(@cols, $desc);
                print $data_json (' ' x 3) . "\"compound\": \"$compound\",\n";
                print $data_json (' ' x 3) . "\"process\": \"$process\",\n";
            }
	}
	
	my $count = $data{$key}->{count};
	my $count_pct = sprintf("%.2f",100* $count/$total_count);
	my $totalAvgDepth = exists $data{$key}->{totalAvgDepth} ? sprintf("%.2f",$data{$key}->{totalAvgDepth}) : 0;
	my $totalAvgDepth_pct = $total_depth != 0 ? sprintf("%.2f",100* $totalAvgDepth/$total_depth) : 0;
	push(@cols, ($count, $count_pct, $totalAvgDepth, $totalAvgDepth_pct));
	
	
	print $data_json (' ' x 3) . "\"count\": $count,\n";
	print $data_json (' ' x 3) . "\"countpct\": $count_pct,\n";
	print $data_json (' ' x 3) . "\"abund\": $totalAvgDepth,\n";
	print $data_json (' ' x 3) . "\"abundpct\": $totalAvgDepth_pct,\n";
	
	my $i = 0; 
	for my $name (sort keys %{$data{$key}}){
	    next if $name eq "count";
	    next if $name eq "totalAvgDepth";
	    my $abund = sprintf("%.2f",$data{$key}->{$name});
	    my $total_abund = $stats{$feature}->{$name};
	    #print STDERR "****************$name, $feature, $abund, $total_abund\n";
	    my $abund_pct = $total_abund > 0 ? sprintf("%.2f",100* $abund/$total_abund) : 0;
	    push(@cols, ($abund, $abund_pct));
	    
#print $data_tab "\t$abund\t$abund_pct";
	    if($i == @msamples - 1){
		print $data_json (' ' x 2) . "\"abund_$name\": $abund,\n";
		print $data_json (' ' x 2) . "\"abundpct_$name\": $abund_pct\n";
	    }
	    else{
		print $data_json (' ' x 2) . "\"abund_$name\": $abund,\n";
		print $data_json (' ' x 2) . "\"abundpct_$name\": $abund_pct,\n";
		
	    }
	    
	    $i++;
	}
	if($index == $keyCount - 1){
	    print $data_json (' ' x 2) . "}\n";
	}
	else{
	    print $data_json (' ' x 2) . "},\n";
	}
	
	print $data_tab join("\t", @cols), "\n";
	$index++;
    }
    
    print $data_json (' ' x 2) . "]\n";
    
    close($data_tab);

    
    


}

sub output_treecol{

    my ($outFile, $firstcol_str, $samples) = @_;
    #msg("Writing to $outFile");
    open (TREECOL,">$outFile")or die "Could not open $outFile to write, $!\n";
    
    my @msamples = @$samples;
     my %treecol = ();
    $treecol{"Count"} = "Count";
    $treecol{"Count_pct"} = "Count%";
    my @colkeys = ("Count","Count_pct");

    if(scalar (@msamples)){
	$treecol{"Abund"} = "Abund";
	$treecol{"Abund_pct"} = "Abund%";
	push(@colkeys, "Abund");
	push(@colkeys, "Abund_pct");
    }
    push(@colkeys, map{"Abund_$_"} @msamples);
    push(@colkeys, map{"Abund_pct_$_"} @msamples);


    foreach my $msample (@msamples){
	$treecol{"Abund_pct_$msample"} = "$msample%";
    }

    print TREECOL "treecol = [\n";
    print TREECOL "{\nheader: \"$firstcol_str\"\n},\n";

    my $size = @colkeys;

    my $i = 0;
    for (; $i < $size; $i++){
	my $key = $colkeys[$i];
	if($i == $size -1){
	    print TREECOL "{width: \"100px\", value: \"$key\", header: \"$treecol{$key}\"}\n" if exists $treecol{$key};
	}
	else{
	    print TREECOL "{value: \"$key\", header: \"$treecol{$key}\"},\n" if exists $treecol{$key};
	}

    }

    print TREECOL "];\n";
    close(TREECOL);

}



sub output_pathway_profile{
    my($sref, $dataref, $feature, $attr, $outdir, $iddesc) = @_;
    my %stats = %$sref;
    my %data = %$dataref;
    my $total_count = $stats{$feature}->{totalCount};
    my $total_depth = exists $stats{$feature}->{totalAvgDepth} ? $stats{$feature}->{totalAvgDepth} : 0;
    
    my @msamples = grep !/totalCount|totalAvgDepth/, sort keys %{$stats{$feature}};
    
    my @header = ("Count", "Count_pct", "Abund", "Abund_pct");
    unshift(@header, "#Pathway_id", "Name", "KOs", "total_family", "total_family_found") if($attr =~ /kegg/);
    unshift(@header, "#Pathway_id", "Type", "Name", "ECs", "total_family", "total_family_found") if($attr =~ /metacyc/);

    foreach (@msamples){
	push(@header, "Abund\_$_");
	push(@header, "Abund\_pct_$_");
    }
    
    open  my $data_tab, ">", "$outdir/$attr.$feature.profile.tab.txt" or die "Could not open $outdir/$attr.$feature.profile.tab.txt to write, $!\n";
    print $data_tab join("\t", @header), "\n";
    open my $data_json, ">", "$outdir/$attr.$feature.profile.datatable.json" or die "Could not open $outdir/$attr.$feature.profile.datatable.json to write, $!\n";

    #output tablecol
    output_tablecols($data_json, \@msamples, $attr);
    
    print $data_json "\ndata = [\n";
    my $keyCount = keys %data;
    
    my $index = 0;
    #my @cols = ();
    for my $key (keys %data){
	my @cols = ();
	push(@cols, $key);
	
	print $data_json (' ' x 2) . "{\n";
	print $data_json (' ' x 3) . "\"pid\": \"$key\",\n";
	
		
	if(defined $iddesc){
	    if($attr eq "kegg"){
		my $name = exists $iddesc->{$key}->{pname} ? $iddesc->{$key}->{pname} : "";
		my $fam_all = exists $iddesc->{$key}->{fam_all} ? $iddesc->{$key}->{fam_all} : "";
		my $fam_found = exists $iddesc->{$key}->{fam_found} ? $iddesc->{$key}->{fam_found} : "";
		my @kos = grep !/pname|fam_found|fam_all/, sort keys %{$keggs->{$key}};
		my $kos_str = join("+", @kos);
		
		push(@cols, ("\"$name\"",$kos_str, $fam_all, $fam_found));
		
		print $data_json (' ' x 3) . "\"name\": \"$name\",\n";
		print $data_json (' ' x 3) . "\"fam_all\": \"$fam_all\",\n";
		print $data_json (' ' x 3) . "\"fam_found\": \"$fam_found\",\n";
		print $data_json (' ' x 3) . "\"kos\": \"$kos_str\",\n";
	    }
	    elsif($attr eq "metacyc"){
		my $name = exists $iddesc->{$key}->{pname} ? $iddesc->{$key}->{pname} : "";
		my $type = exists $iddesc->{$key}->{ptype} ? $iddesc->{$key}->{ptype} : "";
		my $fam_all = exists $iddesc->{$key}->{fam_all} ? $iddesc->{$key}->{fam_all} : "";
		my $fam_found = exists $iddesc->{$key}->{fam_found} ? $iddesc->{$key}->{fam_found} : "";
		
		my @ecs = grep !/pname|fam_found|fam_all|ptype/, sort keys %{$metacycs->{$key}};
		my $ecs_str = join("+", @ecs);
		
		push(@cols, ("\"$type\"", "\"$name\"",$ecs_str, $fam_all, $fam_found));
		print $data_json (' ' x 3) . "\"type\": \"$type\",\n";
		print $data_json (' ' x 3) . "\"name\": \"$name\",\n";
		
		print $data_json (' ' x 3) . "\"fam_all\": \"$fam_all\",\n";
		print $data_json (' ' x 3) . "\"fam_found\": \"$fam_found\",\n";
	    }
	    
	}
	
	my $count = $data{$key}->{count};
	my $count_pct = sprintf("%.2f",100* $count/$total_count);
	my $totalAvgDepth = exists $data{$key}->{totalAvgDepth} ? sprintf("%.2f",$data{$key}->{totalAvgDepth}) : 0;
	my $totalAvgDepth_pct = $total_depth != 0 ? sprintf("%.2f",100* $totalAvgDepth/$total_depth) : 0;
	push(@cols, ($count, $count_pct, $totalAvgDepth, $totalAvgDepth_pct));
	
	
	print $data_json (' ' x 3) . "\"count\": $count,\n";
	print $data_json (' ' x 3) . "\"countpct\": $count_pct,\n";
	print $data_json (' ' x 3) . "\"abund\": $totalAvgDepth,\n";
	print $data_json (' ' x 3) . "\"abundpct\": $totalAvgDepth_pct,\n";
	
	my $i = 0; 
	for my $name (sort keys %{$data{$key}}){
	    next if $name eq "count";
	    next if $name eq "totalAvgDepth";
	    next if $name eq "pname";
	    my $abund = sprintf("%.2f",$data{$key}->{$name});
	    my $total_abund = $stats{$feature}->{$name};
	    my $abund_pct = $total_abund > 0 ? sprintf("%.2f",100* $abund/$total_abund) : 0;
	    #my $abund_pct = sprintf("%.2f",100* $abund/$total_abund);
	    push(@cols, ($abund, $abund_pct));
	    
#print $data_tab "\t$abund\t$abund_pct";
	    if($i == @msamples - 1){
		print $data_json (' ' x 2) . "\"abund_$name\": $abund,\n";
		print $data_json (' ' x 2) . "\"abundpct_$name\": $abund_pct\n";
	    }
	    else{
		print $data_json (' ' x 2) . "\"abund_$name\": $abund,\n";
		print $data_json (' ' x 2) . "\"abundpct_$name\": $abund_pct,\n";
		
	    }
	    
	    $i++;
	}
	if($index == $keyCount - 1){
	    print $data_json (' ' x 2) . "}\n";
	}
	else{
	    print $data_json (' ' x 2) . "},\n";
	}
	
	print $data_tab join("\t", @cols), "\n";
	$index++;
    }
    
    print $data_json (' ' x 2) . "]\n";
    
    close($data_tab);

    
    


}


sub output_tablecols{

    my ($jsonFileHandler, $samples, $attr) = @_;
    my @msamples = @$samples;
    
    print $jsonFileHandler "tablecol = [\n";
    
    if($attr eq "taxon"){
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Taxon\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"taxon\"\n";
	print $jsonFileHandler "},\n";
    }
    elsif($attr eq "kegg"){
        print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"KEGG_ID\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"pid\",\n";
	print $jsonFileHandler (' ' x 2), "\"render\": function(data,type,full,meta){\n";
	print $jsonFileHandler (' ' x 4), "pid = \'ko\' + full.pid;\n";
	print $jsonFileHandler (' ' x 4), "return \'<a href=\"https://www.genome.jp/dbget-bin/www_bget?pathway+\' + pid + \'\">\' + pid + \'</a>\';\n";
	print $jsonFileHandler (' ' x 2), "}\n";

        print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Name\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"name\"\n";
        print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Fam_all\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"fam_all\"\n";
        print $jsonFileHandler "},\n";

	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Fam_found\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"fam_found\",\n";
	print $jsonFileHandler (' ' x 2), "\"render\": function(data,type,full,meta){\n";
	print $jsonFileHandler (' ' x 4), "pid = \'ko\' + full.pid;\n";
	print $jsonFileHandler (' ' x 4), "kos = full.kos;\n";
	print $jsonFileHandler (' ' x 4),"if(data > 0){\n";
	print $jsonFileHandler (' ' x 6), "return \'<a href=\"https://www.genome.jp/kegg-bin/show_pathway?\' + pid + \'+\' + kos + \'\">\' + data + \'</a>\';\n";
	print $jsonFileHandler (' ' x 4),"}\n";
	print $jsonFileHandler (' ' x 4),"else{\n";
	print $jsonFileHandler (' ' x 6),"return data;\n";
	print $jsonFileHandler (' ' x 4),"}\n";
	print $jsonFileHandler (' ' x 2),"}\n";
	
        print $jsonFileHandler "},\n";
    }
    elsif($attr eq "metacyc"){
        print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"MetaCyc_ID\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"pid\",\n";
	print $jsonFileHandler (' ' x 2), "\"render\": function(data,type,full,meta){\n";
	print $jsonFileHandler (' ' x 4), "return \'<a href=\"https://biocyc.org/META/NEW-IMAGE?type=PATHWAY&object=\' + data + \'\">\' + data + \'</a>\';\n";
	print $jsonFileHandler (' ' x 2), "}\n";
	print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Type\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"type\"\n";
        print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Name\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"name\"\n";
        print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Fam_all\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"fam_all\"\n";
        print $jsonFileHandler "},\n";

	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Fam_found\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"fam_found\"\n";
	print $jsonFileHandler "},\n";
    }
    
    elsif($attr eq "metabolic"){
        	
	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Compound\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"compound\"\n";
	print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Process\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"process\"\n";
	print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
        print $jsonFileHandler (' ' x 2), "\"title\": \"Gene\",\n";
        print $jsonFileHandler (' ' x 2), "\"data\": \"gene\"\n";
        print $jsonFileHandler "},\n";

    }
    elsif($attr eq "pfam"){
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Accession\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"accession\"\n";
	print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"ID\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"id\",\n";
	print $jsonFileHandler (' ' x 2), "\"render\": function(data,type,full,meta){\n";
	print $jsonFileHandler (' ' x 4), "return \'<a href=\"https://pfam.xfam.org/family/' + data + \'\">\' + data + \'</a>\';\n";
	print $jsonFileHandler (' ' x 2), "}\n";
	
	print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Description\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"desc\"\n";
	print $jsonFileHandler "},\n";
	
	
    }
    elsif($attr eq "tigrfam"){
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Accession\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"accession\",\n";
	
	print $jsonFileHandler (' ' x 2), "\"render\": function(data,type,full,meta){\n";
	print $jsonFileHandler (' ' x 4), "return \'<a href=\"http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=' + data + \'\">\' + data + \'</a>\';\n";
	print $jsonFileHandler (' ' x 2), "}\n";
	print $jsonFileHandler "},\n";

	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Name\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"name\"\n";
	print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Function\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"desc\"\n";
	print $jsonFileHandler "},\n";
    }
    elsif($attr eq "ko"){
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"KO\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"koid\",\n";

	
	print $jsonFileHandler (' ' x 2), "\"render\": function(data,type,full,meta){\n";
	print $jsonFileHandler (' ' x 4), "return \'<a href=\"https://www.genome.jp/dbget-bin/www_bget?ko\:' + data + \'\">\' + data + \'</a>\';\n";
	print $jsonFileHandler (' ' x 2), "}\n";

	print $jsonFileHandler "},\n";

	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Name\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"name\"\n";
	print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Definition\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"definition\"\n";
	print $jsonFileHandler "},\n";
    }
    elsif($attr eq "ec"){
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"EC\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"ecid\",\n";
	
	print $jsonFileHandler (' ' x 2), "\"render\": function(data,type,full,meta){\n";
	print $jsonFileHandler (' ' x 4), "return \'<a href=\"https://enzyme.expasy.org/EC/' + data + \'\">\' + data + \'</a>\';\n";
	print $jsonFileHandler (' ' x 2), "}\n";
		
	print $jsonFileHandler "},\n";

	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Name\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"definition\"\n";
	print $jsonFileHandler "},\n";
    }
    elsif($attr eq "go"){
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"GO\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"goid\",\n";

	print $jsonFileHandler (' ' x 2), "\"render\": function(data,type,full,meta){\n";
	print $jsonFileHandler (' ' x 4), "return \'<a href=\"http://amigo.geneontology.org/amigo/term/' + data + \'\">\' + data + \'</a>\';\n";
	print $jsonFileHandler (' ' x 2), "}\n";
		
	print $jsonFileHandler "},\n";

	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Name\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"name\"\n";
	print $jsonFileHandler "},\n";
    }
    
    print $jsonFileHandler "{\n";
    print $jsonFileHandler (' ' x 2), "\"title\": \"Count\",\n";
    print $jsonFileHandler (' ' x 2), "\"data\": \"count\"\n";
    print $jsonFileHandler "},\n";

    print $jsonFileHandler "{\n";
    print $jsonFileHandler (' ' x 2), "\"title\": \"Count%\",\n";
    print $jsonFileHandler (' ' x 2), "\"data\": \"countpct\"\n";


    if(@msamples){
	print $jsonFileHandler "},\n";

	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Abund\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"abund\"\n";
	print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Abund%\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"abundpct\"\n";
	print $jsonFileHandler "},\n";
	
	
	
	for my $i (0 .. $#msamples){
	    print $jsonFileHandler "{\n";
	    
	    print $jsonFileHandler (' ' x 2),"\"title\": \"Abund\_$msamples[$i]\",\n";
	    print $jsonFileHandler (' ' x 2),"\"data\": \"abund\_$msamples[$i]\"\n";
	    print $jsonFileHandler "},\n";
	    
	    print $jsonFileHandler "{\n";
	    print $jsonFileHandler (' ' x 2),"\"title\": \"Abund\_$msamples[$i]%\",\n";
	    print $jsonFileHandler (' ' x 2),"\"data\": \"abundpct\_$msamples[$i]\"\n";

	    if($i == $#msamples){
		print $jsonFileHandler "}\n";
	    }
	    else{
		print $jsonFileHandler "},\n";
	    }
	}
    }
    else{
	print $jsonFileHandler "}\n";
    }
    print $jsonFileHandler "];\n";
}
sub output_gff{

    my ($seqHash,$outdir) = @_;
    my $gffver = 3;
    
    msg("Writing master.gff.txt file to $outdir");
    open my $master_gff_fh, '>', "$outdir/master.gff.txt";
    my $master_gff_factory = Bio::Tools::GFF->new(-gff_version=>$gffver);
    print $master_gff_fh "##gff-version $gffver\n";
    
    for my $sid (sort {$seqHash{$b}{DNA}->length <=> $seqHash{$a}{DNA}->length} keys %$seqHash) {
	for my $f ( sort { $a->start <=> $b->start } @{ $seqHash->{$sid}{FEATURE} }) {
	    
	    foreach my  $ftag ($f->all_tags()){
		if($ftag =~ /foam_ecs|foam_kos|foam_target|sprot_kegg|sprot_kos|sprot_ec|sprot_go|pfam_go|tigrfam_ec|tigrfam_mainrole|tigrfam_sub1role|_target/){
		    $f->remove_tag($ftag);
		}
	    }
	    print $master_gff_fh $f->gff_string($master_gff_factory),"\n";
        }
    }
    close($master_gff_fh);
}
