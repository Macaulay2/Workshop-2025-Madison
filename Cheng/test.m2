restart
R = QQ[x,y]
J = ideal(x^2 + y^2, x^2 + y^2 + x*y - 1)
a = x
B := basis(R/J)

R := ring J;
FF := coefficientRing R;
--MO := if not instance(o.MonomialOrder, Nothing) then o.MonomialOrder else (options R).MonomialOrder;
S := newRing(R);
F := sub(J, S);
G := gb(F, ChangeMatrix => true);
aS := sub(a, S);
BS := sub(B, S);
P := last coefficients(BS%F); -- change of basis matrix
V := aS * BS - lift(aS * sub(BS * (inverse P), S/F), S) * P
gens G
HVG := V // gens G;
HGF := getChangeMatrix G;
assert(gens G * HVG - V == 0);
assert(gens F * HGF - gens G == 0);
preH0 := HGF * HVG;
H0 = sub(preH0, ring J);
numcols H0


--H0 := getH0(x, J, Strategy => null)
--H1 := lift(basis(R/J), R);
H1 := sub(syz(gens G), ring J) ;
Theta = random(QQ^(numcols H1), QQ^(numcols H0));
H0 + H1*Theta