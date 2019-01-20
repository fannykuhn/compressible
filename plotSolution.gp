n=1000
reset
set term gif animate delay 10
set output "pression_ordre1.gif"
set cbrange[-0.0:3.0]

do for [j=0:n] {
  plot[-0.01:6.29][-0.01:6.29] sprintf("result%i_p.txt",j) u 1:2:11 with image title sprintf("%i",j)
;}
