use warnings;
use strict;

my $input = shift;
my $input2 = shift;
my %hash;
my $n = 0;
my $e = 0;

open IN,$input or die $!;
while(<IN>){
    chomp;
    $n++;
    if ($n == 1){
#        print "$_\t";
    }
    my @all = split(/\t/,$_);
    my $name = shift@all;
    my $rest = join("\t",@all);
    if(exists $hash{$name}){
        print "something is wrong\n";
    }else{
        $hash{$name} = $rest;
    }
}
close IN;

open IN,$input2 or die $!;
while(<IN>){
    chomp;
    $e++;
    if($e == 1){
#        print "$_\n";
    }
    my @all = split(/\t/,$_);
    my $name = shift@all;
    my $rest = join("\t",@all);
    if(exists $hash{$name}){
         print "$name\t$hash{$name}\t$rest\n";
    }
}
close IN;


