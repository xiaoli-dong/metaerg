
#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Time::Piece;
use HTML::Entities;
use Bio::Tools::GFF;
use FindBin;
use lib "$FindBin::Bin";
require "util.pl";

my ($gff,$outdir, $kegg_profile, $metacyc_profile,$fprot);

$outdir = ".";

&GetOptions(
    "g=s" =>\$gff,
    "k=s" =>\$kegg_profile,
    "m=s" =>\$metacyc_profile,
    "fp=s" =>\$fprot,
    "d=s" =>\$outdir
    );

($gff) ||
    die "Name:\n".
    "  $0 0.1 by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "  Generate gene to function mapping profiles\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -g <gff file>\n".
    "  -k <kegg pathway profile>\n".
    "  -m <metacyc pathway profile>\n".
    "  --fp <CDs proteomics expression file>\n".
    "  -d <output dir>\n";

my $bin = "$FindBin::RealBin";
#my ($msamples, $contig2depth) = parse_depth_file($fdepth);
my %taxon2genes_cds = ();
my %taxon2genes_ssu = ();
my %taxon2genes_lsu = ();
my %pfam2genes = ();
my %pfam = ();
my %tigrfam2genes = ();
my %tigrfam = ();
my %stats = ();
my %metabolic2genes = ();
my %metabolic = ();

my %keggs = ();
my %keggPathway2genes = ();
my %metacycPathway2genes = ();

open(KEGG, $kegg_profile) or die "Could not open $kegg_profile to read, $!\n";
while(<KEGG>){
    chomp;
    next if /^#/;
    my @l = split(/\t/, $_);
    my $pid = $l[0];
    my $pname = $l[1];
    $pname =~ s/^\"|"$//g;
    
    my $fam_all = $l[2];
    my $fam_found = $l[3];
    my $ko = $l[4];
    my $koname = $l[5];
    $koname =~ s/^\"|\"$//g;
    my $hits = $l[6];
    
    $keggs{$pid}->{pname} = $pname;
    $keggs{$pid}->{fam_all} = $fam_all;
    $keggs{$pid}->{fam_found} = $fam_found;
    $keggs{$pid}->{$ko}->{name} = $koname;
    $keggs{$pid}->{$ko}->{hits} = $hits;
    
}


open my $GFF_FILE, "$gff";
my $gff_bio = Bio::Tools::GFF->new(-fh => $GFF_FILE, -gff_version => 3);
while (my $f = $gff_bio->next_feature) {
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
    
    my $mdepth_col_count = $f->get_tag_values("mdepth_cols");
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
close($GFF_FILE);



output_profile(\%stats, \%taxon2genes_cds, "cds", "taxon");
output_profile(\%stats, \%taxon2genes_ssu, "ssu", "taxon");
output_profile(\%stats, \%taxon2genes_lsu, "lsu", "taxon");
output_profile(\%stats, \%pfam2genes, "cds", "pfam", \%pfam);
output_profile(\%stats, \%tigrfam2genes, "cds", "tigrfam", \%tigrfam);
output_profile(\%stats, \%metabolic2genes, "cds", "metabolic", \%metabolic);
output_pathway_profile(\%stats, \%keggPathway2genes, "cds", "kegg", \%keggs);


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




sub output_profile{
    my($sref, $dataref, $feature, $attr, $iddesc) = @_;
    my %stats = %$sref;
    my %data = %$dataref;
    my $total_count = $stats{$feature}->{totalCount};
    my $total_depth = $stats{$feature}->{totalAvgDepth};
    my @msamples = grep !/totalCount|totalAvgDepth/, sort keys %{$stats{$feature}};
    
    my @header = ("Count", "Count_pct", "Abund", "Abund_pct");
    unshift(@header, "#Taxon") if($attr =~ /taxon/);
    unshift(@header, "#Accession", "ID", "Description") if $attr =~ /pfam/;
    unshift(@header, "#Accession", "Name", "Function") if $attr =~ /tigrfam/;
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
	my $totalAvgDepth = sprintf("%.2f",$data{$key}->{totalAvgDepth});
	my $totalAvgDepth_pct = sprintf("%.2f",100* $totalAvgDepth/$total_depth);
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
	    my $abund_pct = sprintf("%.2f",100* $abund/$total_abund);
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
        print $jsonFileHandler (' ' x 2), "\"data\": \"pid\"\n";
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
	print $jsonFileHandler (' ' x 2), "\"data\": \"id\"\n";
	print $jsonFileHandler "},\n";
	
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Description\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"desc\"\n";
	print $jsonFileHandler "},\n";
	
	
    }
    elsif($attr eq "tigrfam"){
	print $jsonFileHandler "{\n";
	print $jsonFileHandler (' ' x 2), "\"title\": \"Accession\",\n";
	print $jsonFileHandler (' ' x 2), "\"data\": \"accession\"\n";
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
sub output_pathway_profile{
    my($sref, $dataref, $feature, $attr, $iddesc) = @_;
    my %stats = %$sref;
    my %data = %$dataref;
    my $total_count = $stats{$feature}->{totalCount};
    my $total_depth = $stats{$feature}->{totalAvgDepth};
    
    my @msamples = grep !/totalCount|totalAvgDepth/, sort keys %{$stats{$feature}};
    
    my @header = ("Count", "Count_pct", "Abund", "Abund_pct");
    unshift(@header, "#Pathway_id", "Name") if($attr =~ /kegg/);
    
    foreach (@msamples){
	push(@header, "Abund\_$_");
	push(@header, "Abund\_pct_$_");
    }
    
    open  my $data_tab, ">", "$outdir/$attr.profile.tab.txt" or die "Could not open $outdir/$attr.profile.tab.txt to write, $!\n";
    print $data_tab join("\t", @header), "\n";
    
    open my $data_json, ">", "$outdir/$attr.profile.datatable.json" or die "Could not open $outdir/$attr.profile.datatable.json to write, $!\n";

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
		my @kos = keys %{$keggs{$key}};
		my $kos_str = join("+", @kos);
		
		push(@cols, ($name, $fam_all, $fam_found, $kos_str));
		
		print $data_json (' ' x 3) . "\"name\": \"$name\",\n";
		print $data_json (' ' x 3) . "\"fam_all\": \"$fam_all\",\n";
		print $data_json (' ' x 3) . "\"fam_found\": \"$fam_found\",\n";
	    }
	    
	}
	
	my $count = $data{$key}->{count};
	my $count_pct = sprintf("%.2f",100* $count/$total_count);
	my $totalAvgDepth = sprintf("%.2f",$data{$key}->{totalAvgDepth});
	my $totalAvgDepth_pct = sprintf("%.2f",100* $totalAvgDepth/$total_depth);
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
	    my $abund_pct = sprintf("%.2f",100* $abund/$total_abund);
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

