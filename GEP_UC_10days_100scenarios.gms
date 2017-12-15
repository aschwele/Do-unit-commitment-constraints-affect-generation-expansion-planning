*GEP.gms
*
$eolcom //
option iterlim=999999999;     // avoid limit on iterations
option reslim=1e6;            // timelimit for solver in sec.
option optcr=0.01;             // gap tolerance
option solprint=ON;           // include solution print in .lst file
option limrow=100;            // limit number of rows in .lst file
option limcol=100;            // limit number of columns in .lst file
//------------------------------------------------------------------------------

//-------------Benders iteration and cut count----------------------------------
Sets
i Benders iteration counter /i1*i30/
count(i) cuts in current Benders master problem
;
Sets
d   days /d1*d10/
g   candidate and existing generating units /g1*g10/
g_e(g) existing generating units /g1*g8/
g_c(g) candidate generating units /g9*g10/
k   candidate and existing wind farms /k1*k2/
k_e(k) existing wind farms /k1/
k_c(k) candidate wind farms /k2/
l   loads   /l1*l4/
h   daily hours /h1*h24/
n   system nodes /n1*n6/
Omega(n,n)   nodes m adjacent to node n
/n1.n2, n2.n1, n1.n3, n3.n1, n2.n3, n3.n2, n2.n4, n4.n2,
n3.n6, n6.n3, n4.n5, n5.n4, n4.n6, n6.n4, n5.n6, n6.n5/
o   generation expansion options  /o1*o5/
s   wind power scenarios /s1*s100/
mapG(g,n) mapping of generators in the system topology
/g1.n1, g2.n1, g3.n4, g4.n5, g5.n2, g6.n2, g7.n3, g8.n3, g9.n2, g10.n6/
mapK(k,n) mapping of wind farms in the system topology /k1.n3, k2.n1/
mapL(l,n) mapping of loads in the system topology /l1.n3, l2.n4, l3.n5, l4.n6/
alias(n,m)
;
Parameters
rho(s) probability of scenario s
sigma(d) weight of representative day d
;
rho(s)=1/card(s);
sigma(d)=365/card(d);
Parameters
V_SH(l) value of load shed for load l /l1 300, l2 300, l3 300, l4 300/
lambda_SU(g) start-up cost of unit g
/g1 175, g2 132, g3 175, g4 107, g5 132, g6 223, g7 283, g8 215, g9 300, g10 100/
C(g) marginal cost of the energy offered by unit g
/g1 11, g2 17, g3 11, g4 23, g5 17, g6 19, g7 19, g8 14, g9 13, g10 25/
P_max(g_e) power capacity of existing unit g
/g1 350, g2 100, g3 76, g4 200, g5 350, g6 197, g7 155, g8 100/
P_min(g_e) minimum power output of existing unit g
/g1 90, g2 10, g3 20, g4 0, g5 30, g6 50, g7 25, g8 50/
R_plus(g_e) maximum ramp-up rate of unit g
/g1 90, g2 20, g3 20, g4 20, g5 70, g6 70, g7 33, g8 65/
R_minus(g_e) maximum ramp-down rate of unit g
/g1 90, g2 20, g3 20, g4 20, g5 70, g6 70, g7 33, g8 65/
R_U(g_e) maximum up-reserve that can be deployed by unit g
/g1 45, g2 10, g3 5, g4 20, g5 30, g6 20, g7 12, g8 15/
R_D(g_e) maximum down-reserve that can be deployed by unit g
/g1 45, g2 10, g3 5, g4 20, g5 30, g6 20, g7 12, g8 15/
I_g(g_c) annualized investment cost of candidate unit g /g9 75000, g10 20000/
I_k(k_c) annualized investment cost of candidate wind farm k /k2 100000/
X_cap(k_e) installed capacity of wind farm k /k1 500/
;

table Load_orig(l,d,h) power consumed by load l in day d in hour h
$ondelim
$include LOADPJM10.csv
$offdelim
;
Parameter
Load(l,d,h);
Load(l,d,h) = Load_orig(l,d,h)*0.05;

table W(k,d,h,s)  wind power production realization of farm k in day d and in hour h and scenario s
$ondelim
$include wind.ten_days.hundred_scenarios.csv
$offdelim
;
Parameter
W_DA_max(k,d,h) maximum power production of wind farm k in day d and in hour h that can be scheduled in day-ahead
;
W_DA_max(k,d,h) = sum(s, rho(s)*W(k,d,h,s))
;

table X_min_opt(g_c,o) minimum power output of option o to build candidate unit g
     o1    o2    o3    o4    o5
g9   0     150   400   600   800
g10  0     50    150   200   300  ;
table X_max_opt(g_c,o) capacity of option o to build candidate unit g
     o1    o2    o3    o4    o5
g9   0     250   500   750   1000
g10  0     200   500   800   1000  ;
table R_plus_opt(g_c,o) maximum ramp-up rate of candidate unit g
     o1    o2    o3    o4    o5
g9   0     50    100   150   200
g10  0     100   250   400   500;
table R_minus_opt(g_c,o) maximum ramp-down rate of candidate unit g
     o1    o2    o3    o4    o5
g9   0     50    100   150   200
g10  0     100   250   400   500;
table P_ini(g,d) initial active power output of unit g in day d
    d1   d2      d3      d4      d5      d6      d7      d8      d9      d10
g1  90   90      90      90      90      90      90      90      90      90
g2
g3
g4
g5
g6  150  150     150     150     150     150     150     150     150     150
g7  155  155     155     155     155     155     155     155     155     155
g8  100  100     100     100     100     100     100     100     100     100
g9
g10 ;
table R_ini(g,d) initial reserve deployed by unit g in day d
   d1    d2      d3      d4      d5      d6      d7      d8      d9      d10
g1 45    45      45      45      45      45      45      45      45      45
g2 0
g3 0
g4
g5
g6
g7
g8 50    50      50      50      50      50      50      50      50      50
g9
g10;
table U_ini(g,d) initial commitment status of unit g in day d
    d1  d2      d3      d4      d5      d6      d7      d8      d9      d10
g1  1    1      1       1       1       1       1       1       1       1
g2
g3
g4
g5
g6  1    1      1       1       1       1       1       1       1       1
g7  1    1      1       1       1       1       1       1       1       1
g8  1    1      1       1       1       1       1       1       1       1
g9
g10;
table F(n,n) capacity of transmission line n to m
     n1          n2      n3      n4      n5      n6
n1   0           500     500     0       0       0
n2   500         0       500     450     0       0
n3   500         500     0       0       0       450
n4   0           450     0       0       500     500
n5   0           0       0       500     0       500
n6   0           0       450     500     500     0;
table X_line(n,n) reactance of transmission line n to m
     n1          n2      n3      n4      n5      n6
n1   0           0.001   0.001   0.001   0.001   0.001
n2   0.001       0       0.001   0.001   0.001   0.001
n3   0.001       0.001   0       0.001   0.001   0.001
n4   0.001       0.001   0.001   0       0.001   0.001
n5   0.001       0.001   0.001   0.001   0       0.001
n6   0.001       0.001   0.001   0.001   0.001   0;

Scalars
RPS Renewable target in percent /0.3/
BigM a large value /1e5/
pi pi /3.1416/
;

//------------------------------------------------------------------------------
//-----------Master problem variables-------------------------------------------
//------------------------------------------------------------------------------
Variables
z_mas cost of generation expansion and day-ahead schedule including expected cost in real-time opeation
alpha(s) expected cost in real-time operation
theta_DA(n,d,h) voltage angle at node n in day d and in hour h at the scheduling stage
Phi1(g_c,d,h)  continuous auxiliary variable
Phi2(g_c,d,h)  continuous auxiliary variable
;
Positive variables
x(k_c) power capacity of candidate wind farm k
x_max(g_c) power capacity of candidate unit g
x_min(g_c) minimum power output of candidate unit g
c_SU(g,d,h) start-up cost of unit g in day d and in hour h
p_DA(g,d,h) power scheduled for unit g in day d and in hour h
w_DA(k,d,h) power scheduled for wind farm k in day d and in hour h
;
Binary variables
u_opt(g_c,o) 1 if option o is selected to build candidate generating unit g
u(g,d,h) 1 if unit g is scheduled to be committed in day d and in hour h
;

alpha.LO(s) = -1e7;
x.LO(k_c)=0;
w_DA.LO(k,d,h)=0;
w_DA.UP(k_e,d,h)=W_DA_max(k_e,d,h)*X_cap(k_e);
theta_DA.UP(n,d,h)=pi;
theta_DA.LO(n,d,h)=-pi;
theta_DA.FX('n1',d,h)=0;
c_SU.LO(g,d,h)=0;

//------------------------------------------------------------------------------
//---------------additional storage parameters----------------------------------
//------------------------------------------------------------------------------
*---------------------additional Parameters for Benders subproblems
Parameter
Load_sub(l,h)
W_sub(k,h)
P_ini_sub(g)
R_ini_sub(g)

min_r(g,d,h)
max_r(g,d,h)
exp_r(g,d,h)
;

//---------------------additional Parameters for storage

 //------------------------Bender storing dual values
Parameters
lambda_x(k_c,d,s,i)
lambda_p_DA(g,d,h,s,i)
lambda_u(g,d,h,s,i)
lambda_Phi1(g_c,d,h,s,i)
lambda_Phi2(g_c,d,h,s,i)
lambda_ramp(g_c,d,s,i)
ubset(i)
lbset(i)
diff(i)

//------------------------Bender storing fixed complicating variables
x_fix(k_c)
p_DA_fix(g,d,h)
u_fix(g,d,h)
Phi1_fix(g_c,d,h)
Phi2_fix(g_c,d,h)

x_fixed(k_c,i)
p_DA_fixed(g,d,h,i)
u_fixed(g,d,h,i)
ramp_fixed(g_c,i)
Phi1_fixed(g_c,d,h,i)
Phi2_fixed(g_c,d,h,i)

//------------------------Bender storing subproblem variables
theta_RT_fixed(n,d,h,s,i)
r_fixed(g,d,h,s,i)
p_SH_fixed(l,d,h,s,i)
w_SP_fixed(k,d,h,s,i)

sub_z(d,s,i)
sub_theta_RT(n,d,h,s)
sub_r(g,d,h,s)
sub_p_SH(l,d,h,s)
sub_w_SP(k,d,h,s)

p_DA_fix_sub(g,h)
u_fix_sub(g,h)
Phi1_fix_sub(g,h)
Phi2_fix_sub(g,h)
ramp_fix_sub(g_c)
;

//------------------------------------------------------------------------------
//-----------Subproblem variables-----------------------------------------------
//------------------------------------------------------------------------------
Free variables
z_sub objective of subproblem (for each day and scenario)
theta_RT_sub(n,h)
r_sub(g,h)
Phi1_sub(g_c,h)
Phi2_sub(g_c,h)
p_DA_sub(g,h)
u_sub(g,h)
ramp_sub(g_c)
x_sub(k_c)
Ramp(g_c)
;
Positive variables
p_SH_sub(l,h)
w_SP_sub(k,h)
;

//------------------------------------------------------------------------------
//-------------subproblem equations---------------------------------------------
//------------------------------------------------------------------------------
Equations
obj_sub,
fix_x,
fix_p_DA,
fix_u,
fix_ramp,
fix_Phi1,
fix_Phi2,
RT_a_sub,
RT_d_sub,
RT_f1_sub, RT_f2_sub,
RT_g1_sub, RT_g2_sub, RT_g11_sub, RT_g22_sub,
RT_h1_sub, RT_h2_sub, RT_h11_sub, RT_h22_sub,
RT_k_sub,
linearRT1_sub, linearRT2_sub
;


obj_sub ..
z_sub =e= ( sum((g,h), C(g)*r_sub(g,h)) + sum((l,h), V_SH(l)*p_SH_sub(l,h))) ;

//fix complicating variables ---------------------------------------------------

fix_x(k_c) .. x_sub(k_c) =e= x_fix(k_c);
fix_p_DA(g,h) .. p_DA_sub(g,h) =e= p_DA_fix_sub(g,h);
fix_u(g,h) .. u_sub(g,h) =e= u_fix_sub(g,h);
fix_ramp(g_c) .. ramp_sub(g_c) =e= ramp_fix_sub(g_c);
fix_Phi1(g_c,h) .. Phi1_sub(g_c,h) =e= Phi1_fix_sub(g_c,h);
fix_Phi2(g_c,h) .. Phi2_sub(g_c,h) =e= Phi2_fix_sub(g_c,h);

//real-time operation ----------------------------------------------------------

RT_a_sub(n,h) .. sum(g$mapG(g,n), p_DA_sub(g,h)+r_sub(g,h))
               + sum(k_e$mapK(k_e,n), W_sub(k_e,h)*X_cap(k_e)-w_SP_sub(k_e,h))
               + sum(k_c$mapK(k_c,n), W_sub(k_c,h)*x_sub(k_c)- w_SP_sub(k_c,h))
               - sum(l$mapL(l,n), load_sub(l,h)-p_SH_sub(l,h)) =e=
      sum(m$(Omega(n,m)), 1/X_line(n,m)*(theta_RT_sub(n,h)-theta_RT_sub(m,h)));

RT_d_sub(k_c,h) .. w_SP_sub(k_c,h) =l= W_sub(k_c,h)*x_sub(k_c);

RT_f1_sub(g_e,h) .. P_min(g_e)*u_sub(g_e,h) =l= r_sub(g_e,h)+p_DA_sub(g_e,h);
RT_f2_sub(g_e,h) .. r_sub(g_e,h)+p_DA_sub(g_e,h) =l= P_max(g_e)*u_sub(g_e,h);

RT_g1_sub(g_e) .. -R_minus(g_e) =l= (p_DA_sub(g_e,'h1')+r_sub(g_e,'h1'))
                                            - (P_ini_sub(g_e)+R_ini_sub(g_e));
RT_g2_sub(g_e) .. (p_DA_sub(g_e,'h1')+r_sub(g_e,'h1'))
                            - (P_ini_sub(g_e)+R_ini_sub(g_e)) =l= R_plus(g_e);

RT_g11_sub(g_c) .. -ramp_sub(g_c) =l= (p_DA_sub(g_c,'h1')+r_sub(g_c,'h1'))
                                            - (P_ini_sub(g_c)+R_ini_sub(g_c));
RT_g22_sub(g_c) .. (p_DA_sub(g_c,'h1')+r_sub(g_c,'h1'))
                           - (P_ini_sub(g_c)+R_ini_sub(g_c)) =l=ramp_sub(g_c);

RT_h1_sub(g_e,h)$(ord(h)>1) .. -R_minus(g_e) =l=
          (p_DA_sub(g_e,h)+r_sub(g_e,h)) - (p_DA_sub(g_e,h-1)+r_sub(g_e,h-1));
RT_h2_sub(g_e,h)$(ord(h)>1) .. (p_DA_sub(g_e,h)+r_sub(g_e,h))
                         - (p_DA_sub(g_e,h-1)+r_sub(g_e,h-1)) =l= R_plus(g_e);

RT_h11_sub(g_c,h)$(ord(h)>1) .. -ramp_sub(g_c) =l=
          (p_DA_sub(g_c,h)+r_sub(g_c,h)) - (p_DA_sub(g_c,h-1)+r_sub(g_c,h-1));
RT_h22_sub(g_c,h)$(ord(h)>1) .. (p_DA_sub(g_c,h)+r_sub(g_c,h))
                       - (p_DA_sub(g_c,h-1)+r_sub(g_c,h-1)) =l= ramp_sub(g_c);

RT_k_sub(n,m,h)$(Omega(n,m)) ..
               1/X_line(n,m)*(theta_RT_sub(n,h)-theta_RT_sub(m,h)) =l= F(n,m);

linearRT1_sub(g_c,h) .. Phi1_sub(g_c,h) =l= r_sub(g_c,h)+p_DA_sub(g_c,h);
linearRT2_sub(g_c,h) .. r_sub(g_c,h)+p_DA_sub(g_c,h) =l= Phi2_sub(g_c,h);


//------------------------------------------------------------------------------
//-----------Master problem equations-------------------------------------------
//------------------------------------------------------------------------------
Equations
obj_mas,
cut,
GE_a, GE_b, GE_c, GE_b1,
GE_e,
DA_a,
DA_c1, DA_c2,
DA_e,
DA_f1, DA_f2, DA_f11, DA_f22,
DA_g1, DA_g2, DA_g11, DA_g22,
DA_j, DA_k, DA_l,
Lin1, Lin2, Lin3, Lin4, Lin5, Lin6, Lin7, Lin8, linearDA1, linearDA2
;

obj_mas ..
z_mas =e= sum(g_c, I_g(g_c)*x_max(g_c))
       +  sum(k_c, I_k(k_c)*x(k_c))
       +  sum(d, sigma(d)* ( sum((g,h), c_SU(g,d,h)+C(g)*p_DA(g,d,h))))
                           + sum(s, rho(s)*alpha(s)) ;

cut(s,count) ..
sum(d, sigma(d)* sub_z(d,s,count)) + (
sum((k_c,d), sigma(d)*lambda_x(k_c,d,s,count) *
                                               (x(k_c) - x_fixed(k_c,count)))+
sum((g,d,h), sigma(d)*lambda_p_DA(g,d,h,s,count) *
                                     (p_DA(g,d,h) - p_DA_fixed(g,d,h,count)))+
sum((g,d,h),  sigma(d)*lambda_u(g,d,h,s,count) *
                                           (u(g,d,h) - u_fixed(g,d,h,count)))+
sum((g_c,d), sigma(d)*lambda_ramp(g_c,d,s,count) *
                                         (Ramp(g_c) - Ramp_fixed(g_c,count)))+
sum((g_c,d,h), sigma(d)*lambda_Phi1(g_c,d,h,s,count) *
                                (Phi1(g_c,d,h) - Phi1_fixed(g_c,d,h,count))) +
sum((g_c,d,h), sigma(d)*lambda_Phi2(g_c,d,h,s,count) *
                                (Phi2(g_c,d,h) - Phi2_fixed(g_c,d,h,count))))
    =l= alpha(s);

//generation expansion ---------------------------------------------------------

GE_a(g_c) .. x_max(g_c) =e= sum(o, u_opt(g_c,o)*X_max_opt(g_c,o));

GE_b(g_c) .. x_min(g_c) =e= sum(o, u_opt(g_c,o)*X_min_opt(g_c,o));

GE_b1(g_c) .. Ramp(g_c) =e= sum(o, u_opt(g_c,o)*R_minus_opt(g_c,o));

GE_c(g_c) .. sum(o, u_opt(g_c,o)) =e= 1;

GE_e .. sum(k_c, x(k_c)) + sum(k_e, X_cap(k_e)) =g=
        RPS*( sum(g_c, x_max(g_c)) + sum(g_e, P_max(g_e)) + sum(k_c, x(k_c))
            + sum(k_e, X_cap(k_e)) );

//first stage/day-ahead --------------------------------------------------------

DA_a(n,d,h) .. sum(g$mapG(g,n), p_DA(g,d,h)) + sum(k$mapK(k,n), w_DA(k,d,h))
             - sum(l$mapL(l,n), load(l,d,h)) =e=
         sum(m$(Omega(n,m)), 1/X_line(n,m)*(theta_DA(n,d,h)-theta_DA(m,d,h)));

DA_c1(g_e,d,h) .. P_min(g_e)*u(g_e,d,h) =l= p_DA(g_e,d,h)
                                               + [min_r(g_e,d,h) *u(g_e,d,h)];
DA_c2(g_e,d,h) .. p_DA(g_e,d,h) + [max_r(g_e,d,h) *u(g_e,d,h)] =l=
                                                        P_max(g_e)*u(g_e,d,h);

DA_e(k_c,d,h) .. w_DA(k_c,d,h) =l= sum(s, rho(s)*W(k_c,d,h,s)*x(k_c));

DA_f1(g_e,d) .. -R_minus(g_e) =l= p_DA(g_e,d,'h1') - P_ini(g_e,d);
DA_f2(g_e,d) .. p_DA(g_e,d,'h1') - P_ini(g_e,d) =l= R_plus(g_e);

DA_f11(g_c,d) .. -Ramp(g_c) =l= p_DA(g_c,d,'h1')- P_ini(g_c,d);
DA_f22(g_c,d) ..                p_DA(g_c,d,'h1') - P_ini(g_c,d) =l= Ramp(g_c);

DA_g1(g_e,d,h)$(ord(h)>1) .. -R_minus(g_e) =l= p_DA(g_e,d,h) - p_DA(g_e,d,h-1);
DA_g2(g_e,d,h)$(ord(h)>1) .. p_DA(g_e,d,h) - p_DA(g_e,d,h-1) =l= R_plus(g_e);

DA_g11(g_c,d,h)$(ord(h)>1) .. -Ramp(g_c) =l= p_DA(g_c,d,h) - p_DA(g_c,d,h-1) ;
DA_g22(g_c,d,h)$(ord(h)>1) ..  p_DA(g_c,d,h) - p_DA(g_c,d,h-1)  =l= Ramp(g_c);

DA_j(n,m,d,h)$(Omega(n,m)) ..
                   1/X_line(n,m)*(theta_DA(n,d,h)-theta_DA(m,d,h)) =l= F(n,m);

DA_k(g,d) .. c_SU(g,d,'h1') =g= lambda_SU(g)*(u(g,d,'h1')-U_ini(g,d));

DA_l(g,d,h)$(ord(h)>1) .. c_SU(g,d,h) =g= lambda_SU(g)*(u(g,d,h)-u(g,d,h-1));

//linearization of equations 2.3b and 2.4e--------------------------------------

Lin1(g_c,d,h) .. -u(g_c,d,h)*BigM =l= Phi1(g_c,d,h);
Lin2(g_c,d,h) ..                      Phi1(g_c,d,h) =l= u(g_c,d,h)*BigM;

Lin3(g_c,d,h) .. -(1-u(g_c,d,h))*BigM =l= (Phi1(g_c,d,h)-x_min(g_c));
Lin4(g_c,d,h) .. (Phi1(g_c,d,h)-x_min(g_c)) =l= (1-u(g_c,d,h))*BigM;

Lin5(g_c,d,h) .. -u(g_c,d,h)*BigM =l= Phi2(g_c,d,h);
Lin6(g_c,d,h) ..                      Phi2(g_c,d,h) =l= u(g_c,d,h)*BigM;

Lin7(g_c,d,h) .. -(1-u(g_c,d,h))*BigM =l= (Phi2(g_c,d,h)-x_max(g_c));
Lin8(g_c,d,h) .. (Phi2(g_c,d,h)-x_max(g_c)) =l= (1-u(g_c,d,h))*BigM;

linearDA1(g_c,d,h) .. Phi1(g_c,d,h) =l= p_DA(g_c,d,h)
                                               + [min_r(g_c,d,h) *u(g_c,d,h)];
linearDA2(g_c,d,h) .. p_DA(g_c,d,h)+ [max_r(g_c,d,h) *u(g_c,d,h)] =l=
                                                                Phi2(g_c,d,h);

//------------------------------------------------------------------------------
//-----------------model definitions--------------------------------------------
//------------------------------------------------------------------------------
model subproblem / obj_sub,
fix_x,
fix_p_DA,
fix_u,
fix_ramp,
fix_Phi1,
fix_phi2,
RT_a_sub,
RT_d_sub,
RT_f1_sub, RT_f2_sub,
RT_g1_sub, RT_g2_sub, RT_g11_sub, RT_g22_sub,
RT_h1_sub, RT_h2_sub, RT_h11_sub, RT_h22_sub,
RT_k_sub,
linearRT1_sub, linearRT2_sub/ ;

model master / obj_mas,  cut ,
GE_a, GE_b, GE_c, GE_b1,
GE_e,
DA_a,
DA_c1, DA_c2,
DA_e,
DA_f1, DA_f2, DA_f11, DA_f22,
DA_g1, DA_g2, DA_g11, DA_g22,
DA_j, DA_k, DA_l,
Lin1, Lin2, Lin3, Lin4, Lin5, Lin6, Lin7, Lin8, linearDA1, linearDA2/ ;

model initial / obj_mas,
GE_a, GE_b, GE_c, GE_b1,
GE_e,
DA_a,
DA_c1, DA_c2,
DA_e,
DA_f1, DA_f2, DA_f11, DA_f22,
DA_g1, DA_g2, DA_g11, DA_g22,
DA_j, DA_k, DA_l,
Lin1, Lin2, Lin3, Lin4, Lin5, Lin6, Lin7, Lin8, linearDA1, linearDA2/;

//------------------------------------------------------------------------------
//------initialization of master problem----------------------------------------
//------------------------------------------------------------------------------

       file opt cplex option file /cplex.opt/;
       put opt;
       put 'lpmethod 6'/;
       put 'startalg 6'/;
       put 'threads 7'/;
       put 'mipemphasis 2'/;
       putclose;

min_r(g,d,h)=  0;
max_r(g,d,h)=  0;
exp_r(g,d,h)=  0;

solve initial minimizing z_mas using MIP;

display
z_mas.l
alpha.l
theta_DA.l
Phi1.l
Phi2.l
x.l
x_max.l
x_min.l
c_SU.l
p_DA.l
w_DA.l
u_opt.l
u.l
;

//------------------------------------------------------------------------------
//-------------------input for Benders algorithm--------------------------------
//------------------------------------------------------------------------------
scalar
ub
lb
tolerance
;
ub = +10e10 ;
lb = -10e10 ;
tolerance = 10 ;

count(i) = no;

//------------------------------------------------------------------------------
//--------------------------initialize matrices as zero-------------------------
//------------------------------------------------------------------------------
Load_sub(l,h)=0;
W_sub(k,h)=0;
P_ini_sub(g)=0;
R_ini_sub(g)=0;

lambda_x(k_c,d,s,i)=0;
lambda_p_DA(g,d,h,s,i)=0;
lambda_u(g,d,h,s,i)=0;
lambda_Phi1(g_c,d,h,s,i)=0;
lambda_Phi2(g_c,d,h,s,i)=0;
lambda_ramp(g_c,d,s,i) =0;

x_fix(k_c)=0;
p_DA_fix(g,d,h)=0;
u_fix(g,d,h) =0;
Phi1_fix(g_c,d,h)=0;
Phi2_fix(g_c,d,h)=0;

x_fixed(k_c,i)=0;
p_DA_fixed(g,d,h,i)=0;
u_fixed(g,d,h,i)=0;
Phi1_fixed(g_c,d,h,i)=0;
Phi2_fixed(g_c,d,h,i)=0;
ramp_fixed(g_c,i)=0;

theta_RT_fixed(n,d,h,s,i)=0;
r_fixed(g,d,h,s,i)=0;
p_SH_fixed(l,d,h,s,i)=0;
w_SP_fixed(k,d,h,s,i)=0;

sub_z(d,s,i)=0;
sub_theta_RT(n,d,h,s)=0;
sub_r(g,d,h,s)=0;
sub_p_SH(l,d,h,s)=0;
sub_w_SP(k,d,h,s)=0;

p_DA_fix_sub(g,h)=0;
u_fix_sub(g,h)=0;
Phi1_fix_sub(g,h)=0;
Phi2_fix_sub(g,h)=0;
ramp_fix_sub(g_c)=0;

ubset(i) =  0;
lbset(i) =  0;
diff(i)  =  10;

//------------------------------------------------------------------------------
//-----------Benders iterations-------------------------------------------------
//------------------------------------------------------------------------------
loop( i $( tolerance > 0.5 ),
count(i) = yes ;

p_DA_fixed(g,d,h,i) = p_DA.L(g,d,h);
u_fixed(g,d,h,i) = u.L(g,d,h);
x_fixed(k_c,i) = x.L(k_c);
Phi1_fixed(g_c,d,h,i) = Phi1.L(g_c,d,h);
Phi2_fixed(g_c,d,h,i) = Phi2.L(g_c,d,h);

p_DA_fix(g,d,h) = p_DA.L(g,d,h);
u_fix(g,d,h) = u.L(g,d,h);
x_fix(k_c) = x.L(k_c);
Phi1_fix(g_c,d,h) = Phi1.L(g_c,d,h);
Phi2_fix(g_c,d,h) = Phi2.L(g_c,d,h);
ramp_fix_sub(g_c) = ramp.L(g_c);

//------------------------------------------------------------------------------
//------solve subproblems-------------------------------------------------------
//------------------------------------------------------------------------------

loop((d,s),

W_sub(k,h) = W(k,d,h,s);
Load_sub(l,h) = Load(l,d,h);
P_ini_sub(g) = P_ini(g,d);
R_ini_sub(g) = R_ini(g,d);

p_DA_fix_sub(g,h) = p_DA_fix(g,d,h);
u_fix_sub(g,h) = u_fix(g,d,h);

Phi1_fix_sub(g_c,h) = Phi1_fix(g_c,d,h);
Phi2_fix_sub(g_c,h) = Phi2_fix(g_c,d,h);

p_SH_sub.LO(l,h)=0;
p_SH_sub.UP(l,h)=Load_sub(l,h);
w_SP_sub.LO(k,h)=0;
w_SP_sub.UP(k_e,h)=W_sub(k_e,h)*X_cap(k_e);
theta_RT_sub.UP(n,h)=pi;
theta_RT_sub.LO(n,h)=-pi;
theta_RT_sub.FX('n1',h)=0;

p_SH_sub.LO(l,h)=0;
p_SH_sub.UP(l,h)=Load_sub(l,h);
w_SP_sub.LO(k,h)=0;
w_SP_sub.UP(k_e,h)=W_sub(k_e,h)*X_cap(k_e);
theta_RT_sub.UP(n,h)=pi;
theta_RT_sub.LO(n,h)=-pi;
theta_RT_sub.FX('n1',h)=0;

solve subproblem using LP minimizing z_sub;

lambda_x(k_c,d,s,i) = fix_x.M(k_c);
lambda_p_DA(g,d,h,s,i) = fix_p_DA.M(g,h);
lambda_u(g,d,h,s,i) = fix_u.M(g,h);
lambda_ramp(g_c,d,s,i) = fix_ramp.m(g_c);
lambda_Phi1(g_c,d,h,s,i) = fix_Phi1.M(g_c,h);
lambda_Phi2(g_c,d,h,s,i) = fix_Phi2.M(g_c,h);

sub_z(d,s,i) = z_sub.L;
sub_theta_RT(n,d,h,s) = theta_RT_sub.L(n,h);
sub_r(g,d,h,s) = r_sub.L(g,h);
sub_p_SH(l,d,h,s) = p_SH_sub.L(l,h);
sub_w_SP(k,d,h,s) = w_SP_sub.L(k,h);

);

min_r(g,d,h) = smin(s,sub_r(g,d,h,s));
max_r(g,d,h) = smax(s,sub_r(g,d,h,s));
exp_r(g,d,h) = sum(s,rho(s)* sub_r(g,d,h,s));

//------------------------------------------------------------------------------
//-----Convergence checking-----------------------------------------------------
//------------------------------------------------------------------------------

ub = sum(g_c, I_g(g_c)*x_max.l(g_c))
   + sum(k_c, I_k(k_c)*x.l(k_c))
   + sum(d, sigma(d)* ( sum((g,h), c_SU.l(g,d,h)+C(g)*p_DA.l(g,d,h))
   + sum(s, rho(s)* ( sum((g,h), C(g)*sub_r(g,d,h,s))
   + sum((l,h), V_SH(l)*sub_p_SH(l,d,h,s))))));

lb = sum(g_c, I_g(g_c)*x_max.l(g_c))
   + sum(k_c, I_k(k_c)*x.l(k_c))
   + sum(d, sigma(d)* ( sum((g,h), c_SU.l(g,d,h)+C(g)*p_DA.l(g,d,h))))
   + sum(s, rho(s)* alpha.l(s));

ubset(i) =  ub;
lbset(i) =  lb;
diff(i)  = {[ub-lb] / ub } * 100 ;
tolerance = abs(((ub-lb)/ub)*100);

if(tolerance < 0.5,
display 'converged';
else

p_DA_fixed(g,d,h,i) = p_DA.L(g,d,h);
u_fixed(g,d,h,i) = u.L(g,d,h);
ramp_fixed(g_c,i) = ramp.L(g_c);
Phi1_fixed(g_c,d,h,i) = Phi1.L(g_c,d,h);
Phi2_fixed(g_c,d,h,i) = Phi2.L(g_c,d,h);
x_fixed(k_c,i) = x.L(k_c);
u_fixed(g,d,h,i) = u.L(g,d,h);

theta_RT_fixed(n,d,h,s,i) = sub_theta_RT(n,d,h,s);
r_fixed(g,d,h,s,i) = sub_r(g,d,h,s);
p_SH_fixed(l,d,h,s,i) = sub_p_SH(l,d,h,s);
w_SP_fixed(k,d,h,s,i) = sub_w_SP(k,d,h,s);

//------------------------------------------------------------------------------
//-----------solve master problem-----------------------------------------------
//------------------------------------------------------------------------------
solve master using MIP minimzing z_mas;
);
);

Parameter
production(g,d,h,s)
windprod(k,d,h,s);

production(g,d,h,s) = p_DA.l(g,d,h) + sub_r(g,d,h,s);
windprod(k_e,d,h,s) =  W(k_e,d,h,s)*X_cap(k_e) - sub_w_SP(k_e,d,h,s);
windprod(k_c,d,h,s) =  W(k_c,d,h,s)*x.l(k_c) - sub_w_SP(k_c,d,h,s);

display
z_mas.L, z_sub.L, alpha.L, ubset, lbset, diff, x.l, u_opt.l, x_max.l, x_min.l,
u.l, c_SU.l, w_DA.l, p_DA.l, theta_DA.l, sub_r, sub_p_SH, sub_w_SP, sub_theta_RT;

execute_unload "10_100.gdx" z_mas.L, z_sub.L, alpha.L, ubset, lbset, diff, x.l,
u_opt.l, x_max.l, x_min.l, u.l, c_SU.l, w_DA.l, p_DA.l, theta_DA.l, sub_r,
sub_p_SH, sub_w_SP, sub_theta_RT, production, windprod
execute 'gdxxrw.exe 10_100.gdx par=z rng=z!'
execute 'gdxxrw.exe 10_100.gdx var=z_mas.l rng=z_mas!'
execute 'gdxxrw.exe 10_100.gdx var=z_sub.l rng=z_sub!'
execute 'gdxxrw.exe 10_100.gdx var=alpha.l rng=alpha!'
execute 'gdxxrw.exe 10_100.gdx par=ubset rng=ubset!'
execute 'gdxxrw.exe 10_100.gdx par=lbset rng=lbset!'
execute 'gdxxrw.exe 10_100.gdx par=diff rng=diff!'
execute 'gdxxrw.exe 10_100.gdx var=x.l rng=x!'
execute 'gdxxrw.exe 10_100.gdx var=u_opt.l rng=u_opt!'
execute 'gdxxrw.exe 10_100.gdx var=x_max.l rng=x_max!'
execute 'gdxxrw.exe 10_100.gdx var=x_min.l rng=x_min!'
execute 'gdxxrw.exe 10_100.gdx var=u.l rng=u!'
execute 'gdxxrw.exe 10_100.gdx var=c_SU.l rng=c_SU!'
execute 'gdxxrw.exe 10_100.gdx var=w_DA.l rng=w_DA!'
execute 'gdxxrw.exe 10_100.gdx var=p_DA.l rng=p_DA!'
execute 'gdxxrw.exe 10_100.gdx var=theta_DA.l rng=theta_DA!'
execute 'gdxxrw.exe 10_100.gdx par=sub_r rng=sub_r!'
execute 'gdxxrw.exe 10_100.gdx par=sub_p_SH rng=sub_p_SH!'
execute 'gdxxrw.exe 10_100.gdx par=sub_w_SP rng=sub_w_SP!'
execute 'gdxxrw.exe 10_100.gdx par=sub_theta_RT rng=sub_theta_RT!'
execute 'gdxxrw.exe 10_100.gdx par=production rng=production!'
execute 'gdxxrw.exe 10_100.gdx par=windprod rng=windprod!'


