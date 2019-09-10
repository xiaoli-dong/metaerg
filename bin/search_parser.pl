#!/usr/bin/env perl

use warnings;
use strict;
use Bio::SeqFeature::Generic;
use Bio::SearchIO;


sub parse_search_results{

    my ($outdir, $seqHash, $coverage, $identity, $hmmevalue) = @_;
    

    msg("Start parsing diamond search against uniprot_sprot db results");
    my $diamond_output = "$outdir/uniprot_sprot.blasttable";
    parse_diamond($diamond_output,"uniprot_sprot", "sprot", $seqHash, $coverage, $identity);
    msg("Finishing parsing diamond search against uniprot_sprot results");


    msg("Start parsing diamond search against genomedb db results");
    $diamond_output = "$outdir/genomedb.blasttable";
    parse_diamond($diamond_output,"genomedb", "genomedb", $seqHash, $coverage, $identity);
    msg("Finish parsing diamond search against genomedb db results");

    msg("Start parsing hmmsearch search against pfam db results");
    my $hmmsearch_output = "$outdir/Pfam-A.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","Pfam-A.hmm", "pfam", $seqHash, $hmmevalue);
    msg("Finish parsing hmmsearch search against pfam db results");

    msg("Start parsing hmmsearch search against tigrfam db results");
    $hmmsearch_output = "$outdir/TIGRFAMs.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","TIGRFAMs.hmm", "tigrfam", $seqHash, $hmmevalue);
    msg("Finish parsing hmmsearch search against tigrfam db results");

    msg("Start parsing hmmsearch search against metabolic_pathway db results");
    $hmmsearch_output = "$outdir/metabolic.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","metabolic.hmm", "homefam", $seqHash, $hmmevalue);
    msg("Finish parsing hmmsearch search against metabolic_pathway db results");

    msg("Start parsing hmmsearch search against casgene db results");
    $hmmsearch_output = "$outdir/casgenes.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","casgenes.hmm", "casgene", $seqHash, $hmmevalue);
    msg("Finish parsing hmmsearch search against casgenes db results");

    msg("Start parsing hmmsearch search against foam hmm db results");
    $hmmsearch_output = "$outdir/FOAM-hmm_rel1a.hmm.hmmer3";
    parse_hmmer3_hmmsearch($hmmsearch_output,"hmmer3","FOAM-hmm_rel1a.hmm", "foam", $seqHash, $hmmevalue);
    msg("Finish parsing hmmsearch search against foam hmm db results");
}

#----------------------------------------------------------------------
sub parse_diamond{

    my ($diamond_output,$dbname, $type, $seqHash, $coverage, $identity) = @_;

#qseqid sseqid qlen qstart qend sstart send qframe pident bitscore evalue qcovhsp
    open(F6, $diamond_output) or die "Could not open $diamond_output to read, !$\n";
    my $num_hit = 0;
    while(<F6>){
	chomp;
	my @line = split(/\t/, $_);
	my ($sid) = $line[0] =~ /^(\S+?)\_cds/;
	for my $f (@{$seqHash->{$sid}{FEATURE}}){
	    if(($f->get_tag_values("ID"))[0] eq $line[0]){
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


sub parse_hmmer3_hmmsearch{
    my ($hmmsearch_out, $format, $dbname, $type, $seqHash, $evalue) = @_;

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

1;

__END__
    
