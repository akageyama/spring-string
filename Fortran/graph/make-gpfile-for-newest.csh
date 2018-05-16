#!/bin/csh -f


set dir = ../src/data
set outfile = ./sgks2.gp

set lastNum = `ls -lt $dir/xyz3.* | head -1 | rev | awk '{print $1}'  | cut -d. -f1 | rev`
echo lastNum = $lastNum

echo \! $0 > $outfile
# echo set xrange \[-10:10\] >> $outfile
# echo set yrange \[-10:10\] >> $outfile
# echo set zrange \[-2:100\] >> $outfile
echo set xrange \[-20:20\] >> $outfile
echo set yrange \[-20:20\] >> $outfile
echo set zrange \[-5:505\] >> $outfile
# echo set view 80, 10 >> $outfile
echo set view 90, 0 >> $outfile

set j = $lastNum

echo splot \'$dir/xyz1.$j\' with linespoints, \'$dir/xyz2.$j\' with linespoints, \'$dir/xyz3.$j\' with linespoints >> $outfile
