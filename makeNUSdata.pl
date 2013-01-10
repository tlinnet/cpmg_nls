#!/usr/bin/perl

if ($#ARGV != 4) {
        print "Delets increments in a 2D fid data file to generate a NUS data set\n";
        print "usage: makeNUSdata.pl <in fid> <NUS table> <out fid> <np> <ni>\n";
        exit;
};

$fidIn = @ARGV[0];
$nusTable = @ARGV[1];
$fidOut = @ARGV[2];
$np = @ARGV[3];
$ni = @ARGV[4];

$bytes = 4*$np;
$i = 0;

open IN, "$nusTable" or die "Cannot open $nusTable for read";
        while(<IN>){;
        @process = split (/\s+/, $_);
        chomp($incrKeep[$i] = @process[0]);
        $i++;
    };
close (IN);

$expectedFileSize = 2*$ni*($bytes+28)+32;
$actualFileSize = -s $fidIn;

if ($expectedFileSize != $actualFileSize){ die "Size of input fid does not match the given np and ni"};


open IN, "<$fidIn" or die "Cannot open $fidIn for read";
open OUT,">$fidOut";


$head="";
$head1="";
$head2="";
$data1="";
$data2="";

sysread(IN,$head,32);
syswrite(OUT,$head,32) ;

for ($i = 0; $i < $ni; $i++){
    sysread(IN,$head1,28);
    sysread(IN,$data1,$bytes);
    sysread(IN,$head2,28);
    sysread(IN,$data2,$bytes);
    if ( $i ~~ @incrKeep){
        syswrite(OUT, $head1, 28);
        syswrite(OUT, $data1, $bytes);
        syswrite(OUT, $head2, 28);
        syswrite(OUT, $data2, $bytes);
    };
};

close(IN);
close(OUT);
