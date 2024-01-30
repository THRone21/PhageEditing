#!/usr/bin/perl
use strict;
my %GC_content; # id=>GC_content
my %sequences; # id=>sequence
my ($id, $base_sum, $gc_sum); # id, 总碱基数，总gc数
my @num; # 中间变量，用于存储单行中某字符的含量

while(my $seq = <>)
{
    chomp($seq);  # 去掉字符串末尾的\n
    if($seq =~ m/^>(.*)/)
    {
        $id = $1;
        next;
    }
    @num = ($seq =~ m/(G|C)/g);
    $GC_content{$id} += @num;
    $gc_sum += @num;    
 
    @num = ($seq =~ m/(.)/g);
    $sequences{$id} += @num;
    $base_sum += @num; 
}


foreach(keys(%GC_content))
{
    if(($GC_content{$_} / $sequences{$_}) > ($gc_sum / $base_sum))
    {
        printf("%s\t%.6f\tbig\n", $_, $GC_content{$_} / $sequences{$_});
    }
    elsif(($GC_content{$_} / $sequences{$_}) < ($gc_sum / $base_sum))
    {
        printf("%s\t%.6f\tsmall\n", $_, $GC_content{$_} / $sequences{$_});
    }
}

