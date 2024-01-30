use strict;
use warnings;

if(!@ARGV){print "<IN1:ID_oldid.txt> <IN2:old_name>\n";exit;} 
open IN1,"$ARGV[0]" or die $!;

my %hash;
while(<IN1>)
{
	next if($.==1);
	chomp;
	my @tmp=split(/\t/,$_,2);
	$hash{$tmp[1]}=$tmp[0];
}
if(exists($hash{$ARGV[1]}))
{
	system("mv $ARGV[1]\.fa $hash{$ARGV[1]}\.fa" );
	}
close IN1;
