
# in this script, sel1 is a selection that is used as reference and sel2 
# must always be within r0 distance from it. If this is not the case, 
# harmonic restraints are used to pull sel2 back to sel1.

# create selections (atom index starts at 1)
# sel1 = miniNLS binding site (heavy atoms)
#set sel1 [addgroup {528 530 532 535 537 538 1056 1058 1060 1063 1064 1066 1068 1070 1072 1074 1075 1112 1114 1116 1119 1120 1122 1124 1125 1126 1128 1130 1132 1134 1135 1165 1167 1169 1171 1173 1177 1178 1212 1214 1216 1220 1221 1222 1224 1226 1229 1231 1232 1233 1235 1238 1239 1240 1242 1244 1246 1248 1252 1253 1254 1256 1258 1261 1263 1264 1297 1299 1301 1303 1305 1309 1310 1666 1668 1670 1673 1676 1677 1678 1681 1682 1709 1711 1713 1716 1717 1719 1721 1722 1723 1725 1727 1729 1731 1732 1769 1771 1773 1776 1777 1778 1781 1782 1819 1821 1823 1826 1827 1828 1829 1830 2363 2365 2367 2370 2371 2372 2375 2376 2410 2412 2414 2417 2418 2420 2422 2423 2424 2426 2428 2430 2432 2433}]
# sel2 = ligand heavy atoms
set sel2 [addgroup {3206 3207 3208 3209 3210 3211 3212 3213 3214 3215 3216 3217 3218 3219 3220 3221 3222 3223 3224 3225 3226 3227 3228 3229 3230 3231 3232 3233 3234 3235 3236}]
#sel3 = peptide heavy atoms
set sel3 [addgroup {3264 3268 3269 3270 3272 3274 3277 3280 3283 3286 3290 3291 3292 3294 3296 3299 3302 3305 3308 3312 3313 3314 3315 3318 3320 3323 3326 3327 3328 3330 3332 3335 3338 3341 3344 3348 3349 3350 3352 3354 3357 3360 3363 3366 3370 3371 3372 3374 3376 3379 3382 3383 3384 3385 3386 3387}]

# spring constant
set k 10

# set target distance r0
set r0 18

# required function
proc calcforces {} {
 # load in atom coordinates (add group computes COM)
 global sel2 k r0 sel3
 loadcoords coord

 # get centers of mass of cm1 and cm2
 #set cm1 [split $coord($sel1) { }]
 set cm1 {-4 8 1.52}
 set cm2 [split $coord($sel2) { }]
 set cm3 [split $coord($sel3) { }]
 #cm1 shifted by 
 
 # extract components of cm1 and cm2
 set x1 [lindex $cm1 0]
 set y1 [lindex $cm1 1]
 set z1 [lindex $cm1 2]

 set x2 [lindex $cm2 0]
 set y2 [lindex $cm2 1]
 set z2 [lindex $cm2 2]

 set x3 [lindex $cm3 0]
 set y3 [lindex $cm3 1]
 set z3 [lindex $cm3 2]

 # compute distance r, lig
 set dx2 [expr $x2 - $x1]
 set dy2 [expr $y2 - $y1]
 set dz2 [expr $z2 - $z1]
 set r2 [expr sqrt($dx2*$dx2 + $dy2*$dy2 + $dz2*$dz2)]

 if {$r2 > $r0} {
 
 # add energy
 addenergy [expr 0.5*$k*($r2-$r0)*($r2-$r0)]

 # add forces
 set fx2 [expr -$k*($r2-$r0)*$dx2/$r2]
 set fy2 [expr -$k*($r2-$r0)*$dy2/$r2]
 set fz2 [expr -$k*($r2-$r0)*$dz2/$r2]
 addforce $sel2 [list $fx2 $fy2 $fz2]
}
 # compute distance r, peptide
 set dx3 [expr $x3 - $x1]
 set dy3 [expr $y3 - $y1]
 set dz3 [expr $z3 - $z1]
 set r3 [expr sqrt($dx3*$dx3 + $dy3*$dy3 + $dz3*$dz3)]
 if {$r3 > $r0} {
 
 # add energy
 addenergy [expr 0.5*$k*($r3-$r0)*($r3-$r0)]

 # add forces
 set fx3 [expr -$k*($r3-$r0)*$dx3/$r3]
 set fy3 [expr -$k*($r3-$r0)*$dy3/$r3]
 set fz3 [expr -$k*($r3-$r0)*$dz3/$r3]
 addforce $sel3 [list $fx3 $fy3 $fz3]
}
}
