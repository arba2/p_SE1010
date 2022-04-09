%% Källa för gemensam grunddata
clear all, close all, clc
SE1010pModellparametrar;   %Läser in SE1010pModellparametrar.m 
FL = (pluft*c*A*v^2)/2;    %Luftmotstånd [N]
%% Lastfall A - Konstant fart v
a = 0;                     %Acceleration
v = v/2;                   %Fart
Hc = 0;
%% Lastfall B - Konstant accelaration a från låg fart
a = a1;                    %Acceleration
v = v/2;                   %Fart
Hc = 0;
%% Lastfall C - Konstant retardation a från maxfart
a = a2;                    %Acceleration
v = v;                     %Fart
Hc = 0;
%% Lastfall D - Kurvtagning
a = 0;                     %Acceleration
v = v/2;                   %Fart
Hc = (m*v^2)/R;            %"Centripetalkraften
%% Gemensam beräknare 
syms Hbi Hby Hf Vbi Vby Vf FD

eq1 = Hc - Hf - Hbi - Hby == 0;                       %x 
eq2 = FD - FL == m*a ;                                %y 
eq3 = Vf + Vbi + Vby - m*g == 0;                      %z
eq4 = FL*h1 + Vf*df + FD*h - (Vbi + Vby)*db == 0;     %Moment kring x i mc
eq5 = (Hf + Hbi + Hby)*h + (Vbi - Vby)*(L/2) == 0;    %Moment kring y i mc
eq6 = Hf*df - (Hbi + Hby)*db == 0;                    %Moment kring z i mc
eq7 = Hbi/Vbi == Hby/Vby;


s = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7], [Hbi, Hby, Hf, Vbi, Vby, Vf, FD]);

Hbi = getfield(s,'Hbi');
Hby = getfield(s,'Hby');
Hf  = getfield(s,'Hf');
Vbi = getfield(s,'Vbi');
Vby = getfield(s,'Vby');
Vf  = getfield(s,'Vf');
FD  = getfield(s,'FD');

%Nedan approximerar bråken som ges av ovan.
Vbi = vpa(Vbi)
Vby = vpa(Vby)
Vf  = vpa(Vf)
Hbi = vpa(Hbi)
Hby = vpa(Hby)
Hf  = vpa(Hf)
FD  = vpa(FD)




%SNITTNING

Mmax_xz_plan = max(abs(Mxz)) %maxvärde
Mmax_xy_plan = max(abs(Mxy)) %maxvärde

M_max = sqrt(Mmax_xz_plan^2 + Mmax_xy_plan^2) %vektoraddera

Tmax_xz_plan = max(abs(Ttot1))
Tmax_xy_plan = max(abs(Ttot2))

T_max = sqrt(Tmax_xz_plan^2 + Tmax_xy_plan^2)


%Stål SS165001 Sträckgräns 640Mpa

Sigma_s = 310*10^6;
n=3;                    %Säkerhetsfaktor
Sigma_eff = Sigma_s/n;

M_tot = M_max;
Mv = -2*FD*dh/2+FK*rd


syms d Wb Wv Sigma_max t_max

qn1 = Wb - (pi*d^3)/32 == 0;
qn2 = Wv - (pi*d^3)/16 == 0;
qn3 = Sigma_max - M_tot/Wb == 0;
qn4 = t_max - Mv/Wv == 0;
qn5 = Sigma_eff - sqrt(abs((Sigma_max)^2+3*(t_max)^2)) == 0;

s = solve([qn1, qn2, qn3, qn4, qn5], [d, Wb, Wv, Sigma_max, t_max])

d = getfield(s,'d')

d = vpa(d)


