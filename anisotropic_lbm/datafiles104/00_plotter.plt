se te jpeg
se pm3d map
se size ratio -1
unset  colorbox
set palette defined (-0.00000 "red", 0.5 "green", 1 "blue")
se cbrange[0:1]
se xrange[0:600]
se yrange[0:600]
unset tics
set cbtics font ", 30"
set cbtics
se ou '00_flwoing_dendrites_20000.jpg'
splot 'phi_20000.dat' us 1:2:3
se te wxt   
