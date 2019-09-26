
#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use HTML::Entities;
use Bio::Tools::GFF;

my ($gff,$feature, $types, $outdir);

$types = "all";
$outdir = ".";
&GetOptions(
    "g=s" =>\$gff,
    "f=s" =>\$feature,
    "t=s" =>\$types,
    "d=s" =>\$outdir
    );

($gff && $feature) ||
    die "Name:\n".
    "  $0 0.1 by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "  Generate gene to function mapping profiles\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -g <gff file>\n".
    "  -t <available annotation type:[all|tRNA|crispr|rRNA|sprot|tigrfam|pfam|casgene|metabolic|genomedb|ec|ko|go], multiple type is separated by colons. For example: sprot:tigrfam:pfam:casgene:metabolic:genomedb:ec:ko:go >\n".
    "  -d <output dir>\n".
    "  -f <feature type, eg.CDS>\n";
    
$types = "tRNA:crispr:rRNA:sprot:tigrfam:pfam:casgene:metabolic:genomedb:ec:ko:go" if $types =~ /all/;

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

	push(@header, "ec_number");
	push(@header, "Taxonomy");
	print {$annot_tab_file_handler{"ec"}} "#", join("\t", (@common_header, @header)), "\n";
    }
    elsif($type eq "ko"){

	push(@header, "ko_number");
	push(@header, "Taxonomy");
	print {$annot_tab_file_handler{"ko"}} "#", join("\t", (@common_header, @header)), "\n";
	
    }
    elsif($type eq "go"){

	push(@header, "go_id");
	push(@header, "Taxonomy");
	print {$annot_tab_file_handler{"go"}} "#", join("\t", (@common_header, @header)), "\n";
	
    }
}


open my $GFF_FILE, "$gff";
my $gff_bio = Bio::Tools::GFF->new(-fh => $GFF_FILE, -gff_version => 3);
while (my $f = $gff_bio->next_feature) {
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
	    my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $ec,$taxon);
	    print {$annot_tab_file_handler{"ec"}} join("\t", @row), "\n";
	    print {$annot_tab_file_handler{"gene2ec"}}  "$featureid\t$ec\n";
	    $ec2genes{$ec}->{$featureid} = 1;
	}
    }
    
    if(exists $annot_tab_file_handler{"ko"} && $f->has_tag("allko_ids")){
	my $count = $f->get_tag_values("allko_ids");
	for (my $i = 0; $i < $count; $i++){
	    my $ko = ($f->get_tag_values("allko_ids"))[$i];
	    my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $ko,$taxon);
	    print {$annot_tab_file_handler{"ko"}} join("\t", @row), "\n";
	    print {$annot_tab_file_handler{"gene2ko"}}  "$featureid\t$ko\n";
	    $ko2genes{$ko}->{$featureid} = 1;
	    
	}
    }
    if(exists $annot_tab_file_handler{"go"} && $f->has_tag("allgo_ids")){
	my $count = $f->get_tag_values("allgo_ids");
	for (my $i = 0; $i < $count; $i++){
	    my $go = ($f->get_tag_values("allgo_ids"))[$i];
	    my @row = ($featureid, $seqid, $start, $end, $len, $strand, $sp, $tm, $go,$taxon);
	    print {$annot_tab_file_handler{"go"}} join("\t", @row), "\n";
	    
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

close($GFF_FILE);
map{close($annot_tab_file_handler{$_})} keys %annot_tab_file_handler;


