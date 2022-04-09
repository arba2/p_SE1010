%% DEL1    Källa för gemensam grunddata
clear all, close all, clc
SE1010pModellparametrar;   %Läser in SE1010pModellparametrar.m 
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
Hc = (m*v^2)/R;            %Centripetalkraften
%% Gemensam beräknare 
FL = (pluft*c*A*v^2)/2;    %Luftmotstånd [N]
syms Hbi Hby Hf Vbi Vby Vf FD       %Sökt

eq11 = Hc - Hf - Hbi - Hby == 0;                       %x 
eq12 = FD - FL == m*a ;                                %y 
eq13 = Vf + Vbi + Vby - m*g == 0;                      %z
eq14 = FL*h1 + Vf*df + FD*h - (Vbi + Vby)*db == 0;     %Moment kring x i mc
eq15 = (Hf + Hbi + Hby)*h + (Vbi - Vby)*(L/2) == 0;    %Moment kring y i mc
eq16 = Hf*df - (Hbi + Hby)*db == 0;                    %Moment kring z i mc
eq17 = Hbi/Vbi == Hby/Vby;


s1 = solve([eq11, eq12, eq13, eq14, eq15, eq16, eq17], [Hbi, Hby, Hf, Vbi, Vby, Vf, FD]);

Hbi = getfield(s1,'Hbi');
Hby = getfield(s1,'Hby');
Hf  = getfield(s1,'Hf');
Vbi = getfield(s1,'Vbi');
Vby = getfield(s1,'Vby');
Vf  = getfield(s1,'Vf');
FD  = getfield(s1,'FD');

%Nedan approximerar bråken som ges av ovan.
Vbi = vpa(Vbi)
Vby = vpa(Vby)
Vf  = vpa(Vf)
Hbi = vpa(Hbi)      %Positiv i negativ riktning
Hby = vpa(Hby)      %Positiv i negativ riktning
Hf  = vpa(Hf)       %Positiv i negativ riktning
FD  = vpa(FD)

%% DEL 2 (endast vänstra kullagret tar upp krafter i sidled)
%Kraftjämvikt
% - pekar i axelns negativa riktning

if FD >= 0               %Fk & Vb enigt momentjämvikt
    Fk = (FD*rh)/rd
    Vb = 0            
else                   
    Fk = 0
    Vb = (FD*rh)/rb
end

Fbi = FD/2              %(FD/(Hbi + Hby))*Hbi   %OM LIVET VORE RÄTTVIST
Fby = Fbi               %(FD/(Hbi + Hby))*Hby

syms Hli Fli Fly Vli Vly

eq21 = Hbi - Hli + Hby == 0;                                                        %Kraftjämvikt i x-led
eq22 = Fbi - Fli + Fk - Fly + Fby == 0;                                             %Kraftjämvikt i y-led
eq23 = Vbi - Vli - Vb - Vly + Vby == 0;                                             %Kraftjämvikt i z-led
%eq24 = - Fbi*rh - Mb*rb + Fk*rd - Fby*rh == 0;                                     %Momentjämvikt kring x
eq25 = Hbi*rh + Vli*b1 + Vb*(L/2 - bb) + Vly*(L - b1) + Hby*rh - Vby*L == 0;        %Momentjämvikt kring y
eq26 = Fbi*0 - Fli*b1 + Fk*(L/2 + bd) - Fly*(L - b1) + Fby*L == 0;                  %Momentjämvikt kring z-axeln

s2 = solve([eq21, eq22, eq23, eq25, eq26], [Hli, Fli, Fly, Vli, Vly]);

Hli = getfield(s2,'Hli');
Fli = getfield(s2,'Fli');
Fly = getfield(s2,'Fly');
Vli = getfield(s2,'Vli');
Vly = getfield(s2,'Vly');


%Nedan approximerar bråken som ges av ovan.
Hli = vpa(Hli)
Fli = vpa(Fli)
Fly = vpa(Fly)
Vli = vpa(Vli)      
Vly = vpa(Vly)      

%%
Snittpunkter = [0, b1, L/2 - bb, L/2 + bd, L - b1, L];

%xz T&M

T = []
M = []
%0<=x<b1
Tz1 = Vbi
My1 = -Hbi*rh - Vbi*x
%b1<=x<(L/2 - bb)
Tz2 = Tz1 - Vli
My2 = My1 + Vli*(x - b1)
%(L/2 -bb)<=x<(L/2 - bd)
Tz3 = Tz2 - Vb
My3 = My2 - Vb*(x - (L/2 - bb))
%(L/2 - bd)<=x<(L - b1)
Tz4 = Tz3
My4 = My3 
%(L - b1)<=x<L
Tz5 = Tz4 - Vly
My5 = My4 - Vly*(x - (L - b1))

%xy T&M
%0<=x<b1
Ty1 = Fbi
Mz1 = Fbi*x
%b1<=x<(L/2 - bb)
Ty2 = Ty1 - Fli
Mz2 = Mz1 + Fli*(x - b1)
%(L/2 -bb)<=x<(L/2 - bd)
Ty3 = Ty2
Mz3 = Mz2 
%(L/2 - bd)<=x<(L - b1)
Ty4 = Ty3 + Fk
Mz4 = Mz3 + Fk*(x - (L/2 - bd))
%(L - b1)<=x<L
Ty5 = Ty4 - Fly
Mz5 = Mz4 - Fly*(x - (L - b1))

%x N&M
%0<=x<b1
Nx1 = - Hbi + Hli - Hby   
Mx1 = FD*rh
%b1<=x<(L/2 - bb)
Nx2 = Nx1 + Hli
Mx2 = Mx1
%(L/2 -bb)<=x<(L/2 - bd)
Nx3 = Nx2
Mx3 = Mx2 - Vb*rb
%(L/2 - bd)<=x<(L - b1)
Nx4 = Nx3
Mx4 = Mx3 + Fk*dr
%(L - b1)<=x<L
Nx5 = Nx4
Mx5 = Mx4 


%% T och m-diagram

Ms = (Hbi*rh) + (Hby*rh)
disp('Inre krafter')
fprintf(' Fk = %g N\n Vli = %g N\n Vly = %g N\n Hli = %g N\n Fli = %g N\n Fly = %g N\n Ms = %g N\n',Fk,Vli,Vly,Hli,Fli,Fly,Ms)

%xy-planet
Txy1 = -FD/2;
Txy2 = -FD/2 + Fli;
Txy3 = Fli - Fk - FD/2;
Txy4 = -FD/2 + Fli - Fk + Fly;

Ttot1 = [Txy1, Txy1, Txy2, Txy2, Txy3, Txy3, Txy4, Txy4]';
figure (1)
intxy = [0, b1, b1, L/2+bd, L/2+bd, L-b1, L-b1, L];
plot(intxy, Ttot1)
title('Tvärkraft XY-planet')
xlabel('Längd [m]')
ylabel('Kraft [N]')

%intervall
xy = [0, b1, b1, L/2+bd, L/2+bd, L-b1, L-b1, L];

My1_1 = (FD/2)*xy(1);
My1_2 = (FD/2)*xy(2);
My2_1 = (FD/2)*xy(3) - Fli*(xy(3)-b1);
My2_2 = (FD/2)*xy(4) - Fli*(xy(4)-b1);
My3_1 = (FD/2)*xy(5) - Fli*(xy(5)-b1) + Fk*(xy(5)-L/2-bd);
My3_2 = (FD/2)*xy(6) - Fli*(xy(6)-b1) + Fk*(xy(6)-L/2-bd);
My4_1 = (FD/2)*xy(7) - Fli*(xy(7)-b1) + Fk*(xy(7)-L/2-bd) - Fly*(xy(7)-(L-b1));
My4_2 = (FD/2)*xy(8) - Fli*(xy(8)-b1) + Fk*(xy(8)-L/2-bd) - Fly*(xy(8)-(L-b1));
Mxy = [My1_1, My1_2, My2_1, My2_2, My3_1, My3_2, My4_1, My4_2];

figure(2)
plot(xy,Mxy)
title('Moment XY-planet')
xlabel('Längd [m]')
ylabel('Moment [Nm]')

%XZ-planet
Txz1 = -Vbi;
Txz2 = -Vbi + Vli;
Txz3 = -Vbi + Vli + Vly;

Ttot2 = [Txz1, Txz1, Txz2, Txz2, Txz3, Txz3]';
figure (3)
intxz = [0, b1, b1, L-b1, L-b1, L];
plot(intxz, Ttot2)
title('Tvärkraft XZ-planet')
xlabel('Längd [m]')
ylabel('Kraft [N]')

xz = [0, b1, b1, L/2+bd, L/2+bd, L-b1, L-b1, L];

Mz1_1 = -Vbi*xz(1) + (Ms/2);
Mz1_2 = -Vbi*xz(2) + (Ms/2);
Mz2_1 = -Vbi*xz(3) + (Ms/2) + Vli*(xz(3)-b1);
Mz2_2 = -Vbi*xz(4) + (Ms/2) + Vli*(xz(4)-b1);
Mz3_1 = -Vbi*xz(5) + (Ms/2) + Vli*(xz(5)-b1);
Mz3_2 = -Vbi*xz(6) + (Ms/2) + Vli*(xz(6)-b1);
Mz4_1 = -Vbi*xz(7) + (Ms/2) + Vli*(xz(7)-b1) + Vly*(xz(7)-L+b1);
Mz4_2 = -Vbi*xz(8) + (Ms/2) + Vli*(xz(8)-b1) + Vly*(xz(8)-L+b1);
Mxz = [Mz1_1, Mz1_2, Mz2_1, Mz2_2, Mz3_1, Mz3_2, Mz4_1, Mz4_2];

figure(4)
plot(xz,Mxz)
title('Moment XZ-planet')
xlabel('Längd [m]')
ylabel('Moment [Nm]')


Mmax_xz_plan = max(abs(Mxz)); %maxvärde
Mmax_xy_plan = max(abs(Mxy)); %maxvärde

M_max = sqrt(Mmax_xz_plan^2 + Mmax_xy_plan^2); %vektoraddera

Tmax_xz_plan = max(abs(Ttot1));
Tmax_xy_plan = max(abs(Ttot2));

T_max = sqrt(Tmax_xz_plan^2 + Tmax_xy_plan^2);


%Stål SS165001 Sträckgräns 640Mpa

Sigma_s = 310*10^6;
n=3;                    %Säkerhetsfaktor
Sigma_eff = Sigma_s/n;

M_tot = M_max;
Mv = -FD*dh+Fk*rd;


syms d Wb Wv Sigma_max t_max

qn1 = Wb - (pi*d^3)/32 == 0;
qn2 = Wv - (pi*d^3)/16 == 0;
qn3 = Sigma_max - M_tot/Wb == 0;
qn4 = t_max - Mv/Wv == 0;
qn5 = Sigma_eff - sqrt(abs((Sigma_max)^2+3*(t_max)^2)) == 0;

s = solve([qn1, qn2, qn3, qn4, qn5], [d, Wb, Wv, Sigma_max, t_max]);

d = getfield(s,'d');

d = vpa(d) 
%Böjskjuvspänningar är ytterst små

Mtot = sqrt(Mxz.^2 + Mxy.^2)'
figure(5)
plot(xz,Mtot)
title('Totala Momentet')
xlabel('Längd [m]')
ylabel('Moment [Nm]')