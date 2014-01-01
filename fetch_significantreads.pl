#!/usr/bin/perl -w
open(RH2,"genecordinatefile_10kb");
#$count=0;
$value1 = 0;
$value2 = 0;
while($line=<RH2>)
{
chomp($line);
@split_line=split(/\t/,$line);
$chr_g=$split_line[0];
$start_g=$split_line[1];
$end_g=$split_line[2];
$gname=$split_line[3];
open(RH1,"CCAT_H3k27me3_H1hesc.significant.region");
while($line1=<RH1>)
{
@split_line1=split(/\t/,$line1);
$chr=$split_line1[0];
$start_p=$split_line1[2];
$end_p=$split_line1[3];
$reads_chip=$split_line1[4];
if(($chr_g eq $chr)&&($start_p<$end_g)&&($end_p>$start_g))
{
    $value1+=$reads_chip;
}
}
print "$gname\t$value1\n";
$value1 = 0;
}
exit;
