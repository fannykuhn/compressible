n=1064
reset
set term gif animate delay 1
set output "choix1_pression_50test.gif"
set cbrange[0.:2.]

do for [j=0:n] {
  plot[-0.01:6.29][-0.01:6.7] sprintf("result%i_p.txt",j) u 1:2:11 with image title sprintf("%i",j)
;}
