#!/bin/bash
# this scripts is to do standard calculation to check
# if the version DMRG-X is right
# jiajunren0522@126.com


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
            cp -r backup "$firstname$name"
            cd "$firstname$name"
            
            if [ $itest -eq 2 ] ; then
                sed -i '/number/i  1' inp
                sed -i '/number/i  1 1 1' inp
                sed -i '/^nstate/c '"3"'' inp
            elif [ $itest -eq 1 ] ; then
                sed -i '/^nstate/c '"1"'' inp
            fi

            sed -i '/^subM/c '"$subM"'' inp
            sed -i '/^spin/c '"$spin"'' inp
            sed -i '/^C2/c '"$C2"'' inp
            sed -i '/^perturbation/c '"$perturbation"'' inp
            sed -i '/^bondorder/c '"$bondorder"'' inp
            sed -i '/^localspin/c '"$localspin"'' inp
            sed -i '/^opmethod/c '"$opmethod"'' inp
            sed -i '/^nstate/c '"$nstate"'' inp
            #qsub DMRG-X.pbs
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
sed -i '/^subM/c '"$subM"'' inp
sed -i '/^spin/c '"$spin"'' inp
sed -i '/^C2/c '"$C2"'' inp
sed -i '/^perturbation/c '"$perturbation"'' inp
sed -i '/^bondorder/c '"$bondorder"'' inp
sed -i '/^localspin/c '"$localspin"'' inp
sed -i '/^opmethod/c '"$opmethod"'' inp
sed -i '/^nstate/c '"$nstate"'' inp
qsub DMRG-X.pbs
cd ..
