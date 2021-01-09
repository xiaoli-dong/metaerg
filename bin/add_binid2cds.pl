#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Tools::GFF;


my ($bindir, $cds, $gff, $prefix);

&GetOptions(
    "d=s" =>\$bindir,
    "c=s" =>\$cds,
    "g=s" =>\$gff,
    "p=s" =>\$prefix
    );

($bindir && $cds && $gff) ||

    die "Name:\n".
    "  $0 0.1 by Xiaoli Dong <xdong\@ucalgary.ca>\n".
    "Synopsis:\n".
    "add bin ids to the protein coding sequence ids in the format of: binid_seqid\n".
    "Usage:\n".
    "  perl $0 \n".
    "  -d <directory name containing all the metabat generated fasta format bin files in the format of \"Bin.*.fa\">\n".
    "  -c <protein coding genes generated from metaerg, which is located in data directory>\n".
    "  -g <master.gff file generated from metaerg, which is located in data directory>\n" .
    "  -p <prefix of the bin file, the preifx will be Bin if your bin files are named as Bin.binid.fa>\n";

$prefix ||= "Bin";

my %contig2bin = ();

opendir(BINDIR, "$bindir") || die "Cannot opendir $bindir: $!\n";
my @dirdata=grep {/$prefix\.\d+\.(fa|fasta)/} readdir(BINDIR);
foreach my $fasta (@dirdata) {
    my ($binid) = $fasta =~ /$prefix\.(\d+)\.(fa|fasta)/;
    print STDERR "found binning file: $fasta\n";
    open (FASTA, "$bindir/$fasta")  || die "Could not open $bindir/$fasta file to read, $!\n";
    $/ = "\n>";
    while(<FASTA>){
	chomp;
	if(my ($seqid,$seq) =  /^>?(\S+).*?\n(.*)/s){
	    #print ">$binid\_$seqid\n$seq\n";
	    $contig2bin{$seqid} = $binid;
	}
    }
    $/ = "\n";
    close(FASTA);
}
close(BINDIR);

my %fid2contig = ();
open my $gff_handle, $gff or die "could not open $gff to read, $!\n";
my $gffio = Bio::Tools::GFF->new(-fh =>$gff_handle, -gff_version => 3);
while (my $f = $gffio->next_feature) {
    my $sid = $f->seq_id;
    my $fid = ($f->get_tag_values("ID"))[0];
    $fid2contig{$fid} = $sid;

}
close($gff_handle);

open (CDS, "$cds") or die "Could not open $cds to read, $!\n";
$/ = "\n>";

while(<CDS>){
    chomp;
    if(/^>?(\S+)(.*?)(OS=.*?)\n(.*)/s){
	my $fid = $1;
	my $fun = $2;
	my $other = $3;
	my $seq = $4;

	$fun =~ s/(^\s+|\s+$)//g;
	$fun =~ s/\s+/_/g;

	my $contigid = $fid2contig{$fid};
	my $binid = "unbinned";
	$binid = $contig2bin{$contigid} if exists $contig2bin{$contigid};

	if($binid eq "unbinned"){
	    print ">$binid\_$fid $fun $other\n";
	}
	else{
	    print ">bin$binid\_$fid $fun $other\n";
	}

	print "$seq\n";
    }
}
$/ = "\n";

close(CDS);
