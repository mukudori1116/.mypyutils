#!/bin/sh

input="com"
prev_tpr="../minimize/min2.gro"
prev_cpt="../minimize/min2.trr"
for (( i=1;i<10; i++));do
    if [ $i -ne 1 ];then
        j=$((i-1))
        prev_tpr=md$j.tpr
        prev_cpt=md$j.cpt
    fi
    if [ ! -f md${i}.gro ];then
        rm -f md$i.out.mdp md$i.tpr
        gmx grompp \
            -f md$i.mdp \
            -c $prev_tpr \
            -r $prev_tpr \
            -p ../top/${input}.top \
            -n ../top/index.ndx \
            -t $prev_cpt \
            -o md$i.tpr \
            -po md$i.out.mdp || exit $?
        rm -f md$i.trr md$i.cpt md$i.gro md$i.edr md$i.log
        gmx mdrun \
            -v -s md$i -deffnm md$i -cpi md$i.cpt
    else
        echo "md${i}.gro present. Start next cycle."
    fi
done
