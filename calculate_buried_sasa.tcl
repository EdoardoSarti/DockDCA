set chA [atomselect top "chain A"]
set chB [atomselect top "chain B"]
set chAB [atomselect top "all"]
set saA [measure sasa 1.4 $chA -samples 1000]
set saB [measure sasa 1.4 $chB -samples 1000]
set saAB [measure sasa 1.4 $chAB -samples 1000]
set BuriedSurfaceArea [expr $saA + $saB - $saAB]
set file [open "tmp_out_sasa.txt" w]
puts $file "$BuriedSurfaceArea"
close $file
exit
