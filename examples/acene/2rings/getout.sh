#!/bin/bash
# this scripts is to do standard calculation to check
# if the version DMRG-X is right
# jiajunren0522@126.com

# single state

touch standard.out

subM="127 127"
perturbation="0 F"
bondorder=1
localspin=1
opmethod="comple"
itest=1
ntest=2
# single state/multistate

for(($itest ; $itest <= $ntest ; itest=$(($itest+1))))
do
    if [ $itest -eq 1 ] ; then
        firstname="single"
        nstate=1
    elif [ $itest -eq 2 ] ; then
        firstname="multi"
        nstate=3
    fi
    echo $firstname
    name=1
    spin=-1
    
    for(($spin ; $spin <= 1 ; spin=$(($spin+1))))
    do
        C2=-1
        for(($C2 ; $C2 <= 1 ; C2=$(($C2+1))))
        do
            cd "$firstname$name"
            
            echo "$firstname$spin$C2" >> ../standard.out
            cat out | grep -A10 "energy converged!" >> ../standard.out
            cat out | grep -A10 "maxiter reached!" >> ../standard.out
            cat out | grep  "RUNTIME" >> ../standard.out
            cat out | grep -A21 "transmoment" >> ../standard.out
            printf '\n\n\n' >> ../standard.out
            
            cd ..
            name=$(($name+1))
        done
    done
done


# Perturbation calculation
spin=0
C2=0
subM="63 127"
perturbation="1 T"
bondorder=0
localspin=0
firstname="perturb"
nstate=1
cp -r backup "$firstname"
cd $firstname
echo "$firstname" >> ../standard.out
cat out | grep -A10 "energy converged!" >> ../standard.out
cat out | grep -A10 "maxiter reached!" >> ../standard.out
cat out | grep  "RUNTIME" >> ../standard.out
printf '\n\n\n' >> ../standard.out
cd ..
