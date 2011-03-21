#!/usr/bin/perl 

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -version => 56
);
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

while($line = <STDIN>){
    chomp($line);
    @t=split("\t", $line);
    # chr	start0	end0	genename	confidence	strand	transcriptname	class
    $chr=$t[0];
    if($chr eq 'chr'){
        print join("\t", @t, "normalized.CpG.content"), "\n";
        next;
    }
    $start=($t[5] eq "-")?$t[1]:($t[1]-1500);
    $end=($t[5] eq "-")?($t[2]+1500):$t[2];

    $chr =~ s/^chr//g;
    
    my $s = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $end, ($t[5] eq '+')?1:-1)->seq;
    
    my @n1 = ($s=~/GC/g); my $n1=@n1;
    my $n2 = length($s);
    my @n3 = (($s=~/C/g) , ($s=~/G/g)); my $n3=@n3;
    
    print join("\t", @t, $s, sprintf("%.3f", ($n1/$n2)/(($n3/$n2/2)*($n3/$n2/2)))), "\n";    
    # Normalized CpG fraction was computed as (observed CpG)/(expected CpG), where expected CpG was calculated as (GC content/2)2. 
}