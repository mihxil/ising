#!/usr/bin/perl
#
# This perl-program can be used to conveniently call the ising program
#
# Michiel Meeuwissen, march 1999


#$seed = 1;
#for($K=0; $K <= 1; $K+=0.1)
#{
#    $q = `ising sx5 sy5 sz5 ntoss100 ni100 i20 seed$seed K$K`;
#    print $q;
#    $seed++;
#}

#test
system("ising sx10 sy10 sa10 mf5 K3.044 seed5 ntoss10000");

