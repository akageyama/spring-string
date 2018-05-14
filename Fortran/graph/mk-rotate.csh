#!/bin/csh -f

set i = 70
set n = 0

while ( $n < 2 ) 
  set j = 0
  while ( $j < 360)
    echo set view $i, $j
    echo replot
    echo pause 0
    @ j+=5
  end
  @ n++
end 
