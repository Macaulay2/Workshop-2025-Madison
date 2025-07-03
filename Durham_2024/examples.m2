uninstallPackage RationalPolytopes

restart
installPackage "RationalPolytopes"


P = convexHull matrix "-1, 2/3"
EP = ehrhartQP P
ES = ehrhartSeries P
value ES
numerator ES
denominator ES

hStarVector P

-- period collapse phenomenon
P = convexHull transpose matrix "1,0; -1,0; 0,1/2; 0,-1/2"
EP = ehrhartQP P
ES = ehrhartSeries P
value ES
numerator ES
denominator ES

hStarVector P
