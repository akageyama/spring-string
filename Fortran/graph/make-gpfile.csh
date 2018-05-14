#!/bin/csh -f


set dir = ../src/data

set lastNum = `ls -lt $dir/xyz3.* | head -1 | rev | awk '{print $1}'  | cut -d. -f1 | rev`

# echo set xrange \[-50:50\]
# echo set yrange \[-50:50\]
# echo set zrange \[-2:100\]
echo set xrange \[-20:20\]
echo set yrange \[-20:20\]
echo set zrange \[-5:505\]
echo set view 80, 10
#echo set view 90, 0

# echo set terminal postscript color

set i = 1
while ( $i < $lastNum )
  if ( $i < 10 ) then
    set j = 000$i
  else if ( $i < 100 ) then
    set j = 00$i
  else
    set j = 0$i
  endif
  
#  echo set output \'$dir/xyz.$j.ps\'
  echo splot \'$dir/xyz1.$j\' with lines, \'$dir/xyz2.$j\' with lines, \'$dir/xyz3.$j\' with lines
  echo pause 0.5
  @ i++
# @ i+=4
end 
