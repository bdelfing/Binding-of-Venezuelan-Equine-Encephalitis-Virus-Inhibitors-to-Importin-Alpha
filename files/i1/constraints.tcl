
# in this script, sel1 is a selection that is used as reference and sel2 
# must always be within r0 distance from it. If this is not the case, 
# harmonic restraints are used to pull sel2 back to sel1.

# create selections (atom index starts at 1)
# sel1 = miniNLS binding site (heavy atoms)
#set sel1 [addgroup {528 530 532 535 537 538 1056 1058 1060 1063 1064 1066 1068 1070 1072 1074 1075 1112 1114 1116 1119 1120 1122 1124 1125 1126 1128 1130 1132 1134 1135 1165 1167 1169 1171 1173 1177 1178 1212 1214 1216 1220 1221 1222 1224 1226 1229 1231 1232 1233 1235 1238 1239 1240 1242 1244 1246 1248 1252 1253 1254 1256 1258 1261 1263 1264 1297 1299 1301 1303 1305 1309 1310 1666 1668 1670 1673 1676 1677 1678 1681 1682 1709 1711 1713 1716 1717 1719 1721 1722 1723 1725 1727 1729 1731 1732 1769 1771 1773 1776 1777 1778 1781 1782 1819 1821 1823 1826 1827 1828 1829 1830 2363 2365 2367 2370 2371 2372 2375 2376 2410 2412 2414 2417 2418 2420 2422 2423 2424 2426 2428 2430 2432 2433}]
# sel2 = ligand heavy atoms
set sel2 [addgroup {3206 3207 3208 3209 3210 3211 3212 3213 3214 3215 3216 3217 3218 3219 3220 3221 3222 3223 3224 3225 3227 3228 3229 3231 3232 3233 3234}]


# spring constant
set k 10

# set target distance r0
set r0 18

# required function
proc calcforces {} {
 # load in atom coordinates (add group computes COM)
 global sel2 k r0
 loadcoords coord

 # get centers of mass of cm1 and cm2
 #set cm1 [split $coord($sel1) { }]
 set cm1 {-4 8 1.52}
 set cm2 [split $coord($sel2) { }]
 #cm1 shifted by 
 
 # extract components of cm1 and cm2
 set x1 [lindex $cm1 0]
 set y1 [lindex $cm1 1]
 set z1 [lindex $cm1 2]

 set x2 [lindex $cm2 0]
 set y2 [lindex $cm2 1]
 set z2 [lindex $cm2 2]

 # compute distance r
 set dx [expr $x2 - $x1]
 set dy [expr $y2 - $y1]
 set dz [expr $z2 - $z1]
 set r [expr sqrt($dx*$dx + $dy*$dy + $dz*$dz)]

 if {$r > $r0} {
 
 # add energy
 addenergy [expr 0.5*$k*($r-$r0)*($r-$r0)]

 # add forces
 set fx [expr -$k*($r-$r0)*$dx/$r]
 set fy [expr -$k*($r-$r0)*$dy/$r]
 set fz [expr -$k*($r-$r0)*$dz/$r]
 addforce $sel2 [list $fx $fy $fz]
}
}
