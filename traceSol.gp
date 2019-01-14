# Output dans un fichier pdf
set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
    size 10,6 # inches

set palette rgbformulae 7,5,15

set view 60,45
#set hidden3d
#set zrange[0:1]
set ztics
unset key

set output "result0.0_p.pdf"
plot "result0.000000_p.txt" u 1:2:3 with image