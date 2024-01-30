#!/usr/bin/perl
use strict;
my %G_content; # id=>G_content
my %C_content; # id=>C_content
my $id; # id, 总碱基数，总gc数
my @num; # 中间变量，用于存储单行中某字符的含量

while(my $seq = <>)
{
    chomp($seq);  # 去掉字符串末尾的\n
    if($seq =~ m/^>(.*)/)
    {
        $id = $1;
        next;
    }
    @num = ($seq =~ m/(G)/g);
    $G_content{$id} += @num;
    @num = ($seq =~ m/(C)/g);
    $C_content{$id} += @num;    
}


foreach(keys(%G_content))
{
    printf("%s\t%.6f\n", $_, ($G_content{$_} - $C_content{$_}) / ($G_content{$_} + $C_content{$_}));
}
