# Output dans un fichier pdf
set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
    size 10,6 # inches

set palette defined (0 0.9 0.9 1, 35 0.3 0.3 1, 50 0.6 0.15 0.4, 70 'red', 100 'yellow')
set cbrange[0.:3.5]
set view 60,45
#set hidden3d
#set zrange[0:1]
set ztics
unset key

set output "ordre1_pression_Bnul_aprem.pdf"
plot "result1776_p.txt" u 1:2:3 with image
