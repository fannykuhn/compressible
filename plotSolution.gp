n=1800
reset
set term gif animate delay 1
set output "nan_MUSCL_choix1.gif"
set cbrange[0.:4.5]
set palette defined (0 0.9 0.9 1, 35 0.3 0.3 1, 50 0.6 0.15 0.4, 70 'red', 100 'yellow')


do for [j=0:n] {
  plot[-0.01:6.29][-0.01:6.7] sprintf("resultMUSCL%i_p.txt",j*10) u 1:2:3 with image title sprintf("%i",j*10)
;}
