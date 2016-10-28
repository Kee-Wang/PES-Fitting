#!/usr/bin/perl -w
#============================================================
# this perl script is used to extract the geometry
# and energy information from the molpro output.
# usually we can extract the info from molpro punch
# file which is much easy to handle. the output is
# a xyz file "mp.mpo" which contains the energy info and
# geometry info
# ============================================================
use strict;
# =========================================================================
# a typical punch file for sp energy calculation mpsp.pun
#=========================================================================#
# ***       SP AND GRAD FOR H5+                                        	  #
# TIME        6-Jan-0  12:30:02 				       	  #
# ATOM  1  H      1.00        -1.97929588    -0.63790518    -0.38329212	  #
# ATOM  2  H      1.00         1.98587069     0.38642153    -0.63450918	  #
# ATOM  3  H      1.00         1.98012288    -0.38590336     0.63396586	  #
# ATOM  4  H      1.00        -1.98333795     0.63860185     0.38538295	  #
# ATOM  5  H      1.00        -0.00335976    -0.00121483    -0.00154752	  #
# RHF STATE 1.1 ENERGY        -2.43918222			       	  #
# RHF STATE 1.1 DIPOLM        -0.00850856     0.00103766    -0.00031234	  #
# MP2 STATE  1.1 ENERGY       -2.51171579			       	  #
# MP2 STATE  1.1 ENERGY       -2.51171579			       	  #
# ---								       	  #
#=========================================================================#

open(SPF,"mpsp.pun") or die "Can not open file mpsp.pun\n";
open(MPO,">mp.mpo")  or die "Can not create file mp.mpo\n";

my @xyz=();
my @symb=();
my @enrg=();

<SPF>;
<SPF>;
while(($_=<SPF>) =~ /ATOM/){
    chomp;
    my @entry=split;
    push(@symb,$entry[2]);
    push(@xyz,@entry[4..6]);
}
my $dim=$#xyz+1;
my $natm=$dim/3;

if ($_ =~ /Energy/){
    chomp;
    my @entry=split;
    push(@enrg,$entry[4]);
}

while(<SPF>){
    next unless /Energy/;
    chomp;
    my @entry=split;
    push(@enrg,$entry[4]);
}
close(SPF);

my $energy=pop(@enrg);

printf MPO "%2d\n",$natm;
printf MPO "%15.8f\n",$energy;

for(my $i=1;$i<=$natm;$i++){
    printf MPO "%2s %15.8f %15.8f %15.8f\n",$symb[$i-1],@xyz[3*$i-3..3*$i];
}

close(MPO);
