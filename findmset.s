#!/bin/tcsh

rm msetresult.dat

set nv = 27

foreach m11 (2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21)
foreach m12 (3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
foreach m13 (4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23)
foreach m14 (5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
foreach m15 (6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25)
foreach m16 (7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26)

if (($m11<$m12)&&($m12<$m13)&&($m13<$m14)&&($m14<$m15)&&($m15<$m16)) then
   echo  "  1   $m11   $m12   $m13   $m14   $m15   $m16   $nv  "  >> msetresult.dat
endif
end
end
end
end
end
end

echo "End" >> msetresult.dat
