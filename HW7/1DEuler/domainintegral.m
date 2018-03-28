clear all
syms ua1 ua2 ub1 ub2 uc1 uc2
syms x g 

ua=ua1+x*ua2;
ub=ub1+x*ub2;
uc=uc1+x*uc2;

p=(g-1)*(uc - ub^2/(2*ua))

F21= int(ub, x)
F22= int(ub^2/ua + p,x)
F23= int(ub*uc/ua  + (ub/ua)*p, x)

F21sol=ccode(simplify(subs(F21,x,1)-subs(F21,x,-1)))
F22sol=ccode(simplify(subs(F22,x,1)-subs(F22,x,-1)))
F23sol=ccode(simplify(subs(F23,x,1)-subs(F23,x,-1)))
