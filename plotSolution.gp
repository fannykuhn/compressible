n=2000
reset
set term gif animate delay 10
set output "choix1_pression_ordre1.gif"
set cbrange[0.:4.5]

do for [j=0:n] {
  plot[-0.01:6.29][-0.01:6.7] sprintf("result%i_p.txt",j*10) u 1:2:11 with image title sprintf("%i",j*10)
;}
