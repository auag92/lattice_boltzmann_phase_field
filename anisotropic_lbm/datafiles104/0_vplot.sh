#!/bin/sh
arg1=$1
>Plot.txt
echo "set term gif animate"	>>Plot.txt
#echo "set view equal xy" >>Plot.txt
echo "set o \"0_velocity.gif\""	>>Plot.txt
echo "unset key"		>>Plot.txt

for i in $(seq 50000 2000 290000)
do
        newname1=$(printf '%s_%s.dat' $arg1 $i)
	echo "plot '$newname1' us 1:2:3:4 w vec" >>Plot.txt
done

