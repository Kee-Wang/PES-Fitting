#!/usr/bin/perl -w
# ============================================================
# this perl script is used to extract the geometry, energy
# and gradient information from the molpro output file.
# usually it is from molpro punch file. the output is a
# xyz file "mp.mpo" which contains the energy,geometry and
# gradient.
# ============================================================
use strict;
#=================================================================================#
# a typical punch file of molpro performing both sp and force                     #
#=================================================================================#
# ***       SP AND GRAD FOR H5+                                                   #
# TIME       31-Dec-0  14:55:48 						  #
# ATOM  1  H      1.00        -1.97929588    -0.63790518    -0.38329212		  #
# ATOM  2  H      1.00         1.98587069     0.38642153    -0.63450918		  #
# ATOM  3  H      1.00         1.98012288    -0.38590336     0.63396586		  #
# ATOM  4  H      1.00        -1.98333795     0.63860185     0.38538295		  #
# ATOM  5  H      1.00        -0.00335976    -0.00121483    -0.00154752		  #
# RHF STATE 1.1 ENERGY        -2.43918222					  #
# RHF STATE 1.1 DIPOLM        -0.00850855     0.00103767    -0.00031234		  #
# MP2 STATE  1.1 ENERGY       -2.51171580					  #
# MP2 STATE  1.1 ENERGY       -2.51171580					  #
# 										  #
# MP2 GRADIENT, ATOM  1:          0.000791299078  -0.003579694473  -0.002135011862#
# MP2 GRADIENT, ATOM  2:         -0.000490439928   0.001683501386  -0.002718730219#
# MP2 GRADIENT, ATOM  3:         -0.000868878710  -0.001594861273   0.002651116643#
# MP2 GRADIENT, ATOM  4:          0.000397342308   0.003719328401   0.002252042261#
# MP2 GRADIENT, ATOM  5:          0.000170677252  -0.000228274041  -0.000049416823#
# ---										  #
#=================================================================================#

open(EGF,"mpeg.pun") or die "Can not open file mpeg.pun\n";
open(MPO,">mp.mpo")  or die "Can not create file mp.mpo\n";

my @xyz=();
my @grad=();
my @symb=();
my @enrg=();

<EGF>;
<EGF>;
while(($_=<EGF>) =~ /ATOM/){
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

until(($_=<EGF>)=~/^$/){
    next unless /Energy/;
    chomp;
    my @entry=split;
    push(@enrg,$entry[4]);
}

while($_=<EGF>){
    if ($_ =~ /GRADIENT/){
        chomp;
        my @entry=split;
        push(@grad,@entry[4..6]);
    }
}
close(EGF);

my $energy=pop(@enrg);

printf MPO "%2d\n",$natm;
printf MPO "%15.8f\n",$energy;

for(my $i=1;$i<=$natm;$i++){
    printf MPO "%2s %15.8f %15.8f %15.8f %19.12f %19.12f %19.12f\n", 
    $symb[$i-1],@xyz[3*$i-3..3*$i-1],@grad[3*$i-3..3*$i-1];
}

close(MPO);
