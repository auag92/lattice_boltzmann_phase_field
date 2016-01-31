#!/bin/sh
arg1=$1
>Plot.txt
echo "set term gif animate"	>>Plot.txt
#echo "set term postscript"	>>phi1.txt
echo "set view map"	>>Plot.txt
echo "set pm3d"		>>Plot.txt
echo "unset surface"	>>Plot.txt
echo "set view 0,0"    >>Plot.txt
echo "set view equal xy" >>Plot.txt
#echo "set pm3d map"	>>ConcA2.txt
echo "set cbrange[0.0:1]" >>Plot.txt
echo "set palette defined (-0.00000 \"red\", 0.5 \"green\", 1 \"blue\")" >>Plot.txt
echo "set o \"0_phi.gif\""	>>Plot.txt
echo "unset key"		>>Plot.txt

#for((int i=6000;i<7000;i=i+1))
for i in $(seq 50000 2000 290000)
do
        newname1=$(printf '%s_%s.dat' $arg1 $i)
	echo "splot '$newname1' u 1:2:3" >>Plot.txt
done

