#!/usr/bin/perl -w
open(RH2,"Genecordinatefile");
$value=0;
$norm=0;
while($line=<RH2>)
{
chomp($line);
@split_line=split(/\t/,$line);
$chr=$split_line[1];
$start_g=$split_line[2];
$end_g=$split_line[3];
$genelength=$end_g-$start_g;
$start=$start_g-1000000;
$end=$end_g+1000000;
open(RH1,"HCNE_locationfile_foreachspecie");
while($line1=<RH1>)
{
@split_line1=split(/\t/,$line1);
$chr1=$split_line1[0];
$start1=$split_line1[1];
$end1=$split_line1[2];
$val=$end1-$start1;

if(($chr1 eq $chr)&&($start1>$start)&&($end1<$end))
{

$value+=$val;
$norm=($value/$genelength)*100;

}
}
printf "%4.4f",$norm;
print "\n";
$value=0;
$norm=0;
}
exit;
