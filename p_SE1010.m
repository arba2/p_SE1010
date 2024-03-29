%% Lastfall A - Konstant fart v
clear all, close all, clc, SE1010pModellparametrar; LF = 'A';
a = 0;                     %Acceleration
v = v/2;                   %Fart
%% Lastfall B - Konstant accelaration a fr�n l�g fart
clear all, close all, clc, SE1010pModellparametrar; LF = 'B';
a = a1;                    %Acceleration
v = v/2;                   %Fart
%% Lastfall C - Konstant retardation a fr�n maxfart
clear all, close all, clc, SE1010pModellparametrar; LF = 'C';
a = a2;                    %Acceleration
v = v;                     %Fart
%% Lastfall D - Kurvtagning
clear all, close all, clc, SE1010pModellparametrar; LF = 'D';
a = 0;                     %Acceleration
v = v/2;                   %Fart
Hc = (m*v^2)/R;            %Centripetalkraften
%% Gemensam j�mviktsber�knare 
FL = (pluft*c*A*v^2)/2;                                %Luftmotst�nd [N]
syms Hbi Hby Hf Vbi Vby Vf FD                          %S�kt
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
%Nedan approximerar br�ken som ges av ovan.
Vbi = vpa(Vbi);
Vby = vpa(Vby);
Vf  = vpa(Vf);
Hbi = vpa(Hbi);      %Positiv i negativ riktning
Hby = vpa(Hby);      %Positiv i negativ riktning
Hf  = vpa(Hf);       %Positiv i negativ riktning
FD  = vpa(FD);
disp('Yttre belastningar')
fprintf(' Vbi = %g N\n Vby = %g N\n Vf = %g N\n Hbi = %g N\n Hby = %g N\n Hf = %g N\n FD = %g N\n ',Vbi,Vby,Vf,Hbi,Hby,Hf,FD)

% DEL 2 %Kraftj�mvikt (endast v�nstra kullagret tar upp krafter i sidled)

% - pekar i axelns negativa riktning

if  FD >= 0              %Fk & Vb enigt momentj�mvikt
    Fk = (FD*rh)/rd;
    Vb = 0;            
else                   
    Fk = 0;
    Vb = (FD*rh)/rb;
end
Fbi = FD/2;              %(FD/(Hbi + Hby))*Hbi   %OM LIVET VORE R�TTVIST
Fby = Fbi;               %(FD/(Hbi + Hby))*Hby

syms Hli Fli Fly Vli Vly
eq21 = Hbi - Hli + Hby == 0;                                                        %Kraftj�mvikt i x-led
eq22 = Fbi - Fli + Fk - Fly + Fby == 0;                                             %Kraftj�mvikt i y-led
eq23 = Vbi - Vli - Vb - Vly + Vby == 0;                                             %Kraftj�mvikt i z-led
%eq24 = - Fbi*rh - Mb*rb + Fk*rd - Fby*rh == 0;   �VERFL�DIG                        %Momentj�mvikt kring x
eq25 = Hbi*rh + Vli*b1 + Vb*(L/2 - bb) + Vly*(L - b1) + Hby*rh - Vby*L == 0;        %Momentj�mvikt kring y
eq26 = Fbi*0 - Fli*b1 + Fk*(L/2 + bd) - Fly*(L - b1) + Fby*L == 0;                  %Momentj�mvikt kring z-axeln

s2 = solve([eq21, eq22, eq23, eq25, eq26], [Hli, Fli, Fly, Vli, Vly]);
Hli = getfield(s2,'Hli');
Fli = getfield(s2,'Fli');
Fly = getfield(s2,'Fly');
Vli = getfield(s2,'Vli');
Vly = getfield(s2,'Vly');
%Nedan approximerar br�ken som ges av ovan.
Hli = vpa(Hli);
Fli = vpa(Fli);
Fly = vpa(Fly);
Vli = vpa(Vli);
Vly = vpa(Vly);      
disp('Inre krafter')
fprintf(' Fk = %g N\n Vb = %g N\n Vli = %g N\n Vly = %g N\n Hli = %g N\n Fli = %g N\n Fly = %g N\n',Fk,Vb,Vli,Vly,Hli,Fli,Fly)


%% SNITTNING ****************************************************************
ix =  [0, b1, b1, L/2-bb, L/2-bb, L/2+bd, L/2+bd, L-b1, L-b1, L]; %Intervallvektor

%Normalkraft i x-riktning
Nx1 = Hbi;
Nx2 = Nx1 - Hli;
Nx3 = Nx2; 
Nx4 = Nx3;
Nx5 = Nx4 + Hby;

Nx = [Nx1, Nx1, Nx2, Nx2, Nx3, Nx3, Nx4, Nx4, Nx5, Nx5]';
figure (1)
plot(ix, Nx)
title(['Normalkraft i x-riktning Lastfall ' [LF]])
xlabel('L�ngd [m]')
ylabel('Kraft [N]')

% Tv�rkraft i y-riktning
Ty1 = -Fbi;
Ty2 = -Fbi + Fli;
Ty3 =  Ty2;
Ty4 = -Fbi + Fli - Fk;
Ty5 = -Fbi + Fli - Fk + Fly;

Ty = [Ty1, Ty1, Ty2, Ty2, Ty3, Ty3, Ty4, Ty4, Ty5, Ty5]';
figure (2)
plot(ix, Ty)
title(['Tv�rkraft i y-riktning Lastfall ' [LF]])
xlabel('L�ngd [m]')
ylabel('Kraft [N]')

% Tv�rkraft i z-riktning
Tz1 = -Vbi;
Tz2 = -Vbi + Vli;
Tz3 = -Vbi + Vli + Vb;
Tz4 =  Tz3;
Tz5 = -Vbi + Vli + Vb + Vly;

Tz = [Tz1, Tz1, Tz2, Tz2, Tz3, Tz3, Tz4, Tz4, Tz5, Tz5]';
figure (3)
plot(ix, Tz)
title(['Tv�rkraft i z-riktning Lastfall ' [LF]])
xlabel('L�ngd [m]')
ylabel('Kraft [N]')

% Moment kring x
Mx1 = -Fbi*rh;
Mx2 = Mx1;
Mx3 = Mx1 - Vb*rb; 
Mx4 = Mx3 + Fk*rd;
Mx5 = Mx4;

Mx = [Mx1, Mx1, Mx2, Mx2, Mx3, Mx3, Mx4, Mx4, Mx5, Mx5]';
figure (4)
plot(ix, Mx)
title(['Moment kring x Lastfall ' [LF]])
xlabel('L�ngd [m]')
ylabel('Kraft [N]')

% Moment kring y
My11 = Hbi*rh - Vbi*ix(1);
My12 = Hbi*rh - Vbi*ix(2);
My21 = Hbi*rh - Vbi*ix(2) + Vli*(ix(2)-ix(2));
My22 = Hbi*rh - Vbi*ix(4) + Vli*(ix(4)-ix(2));
My31 = Hbi*rh - Vbi*ix(4) + Vli*(ix(4)-ix(2)) + Vb*(ix(4)-ix(4));
My32 = Hbi*rh - Vbi*ix(6) + Vli*(ix(6)-ix(2)) + Vb*(ix(6)-ix(4));
My41 = Hbi*rh - Vbi*ix(6) + Vli*(ix(6)-ix(2)) + Vb*(ix(6)-ix(4));  % My4 INTE ETT SNITT 
My42 = Hbi*rh - Vbi*ix(8)  + Vli*(ix(8) -ix(2)) + Vb*(ix(8) -ix(4)) + Vly*(ix(8) -ix(8));
My51 = Hbi*rh - Vbi*ix(8)  + Vli*(ix(8) -ix(2)) + Vb*(ix(8) -ix(4)) + Vly*(ix(8) -ix(8));
My52 = Hbi*rh - Vbi*ix(10) + Vli*(ix(10)-ix(2)) + Vb*(ix(10)-ix(4)) + Vly*(ix(10)-ix(8));

My = [My11, My12, My21, My22, My31, My32, My41, My42, My51, My52];
%ixMy = ix; ixMy(7:8) = [];                  %Kortar ner ix till 4 intervall
figure(5)
plot(ix,My)
title(['Moment kring y Lastfall ' [LF]])
xlabel('L�ngd [m]')
ylabel('Moment [Nm]')

% Moment kring z-axeln
Mz11 = Fbi*ix(1);
Mz12 = Fbi*ix(2);
Mz21 = Fbi*ix(2)  - Fli*(ix(2) -ix(2));
Mz22 = Fbi*ix(4)  - Fli*(ix(4) -ix(2)); 
Mz31 = Fbi*ix(4)  - Fli*(ix(4) -ix(2));  % Mz3 INTE ETT SNITT 
Mz32 = Fbi*ix(6)  - Fli*(ix(6) -ix(2)) + Fk*(ix(6) -ix(6));
Mz41 = Fbi*ix(6)  - Fli*(ix(6) -ix(2)) + Fk*(ix(6) -ix(6));
Mz42 = Fbi*ix(8)  - Fli*(ix(8) -ix(2)) + Fk*(ix(8) -ix(6));
Mz51 = Fbi*ix(8)  - Fli*(ix(8) -ix(2)) + Fk*(ix(8) -ix(6)) - Fly*(ix(8) -ix(8));
Mz52 = Fbi*ix(10) - Fli*(ix(10)-ix(2)) + Fk*(ix(10)-ix(6)) - Fly*(ix(10)-ix(8));

Mz   = [Mz11, Mz12, Mz21, Mz22, Mz31, Mz32 Mz41, Mz42, Mz51, Mz52];
%ixMz = ix; ixMz(5:6) = [];                  %Kortar ner ix till 4 intervall
figure(6)
plot(ix,Mz)             
title(['Moment kring z Lastfall ' [LF]])
xlabel('L�ngd [m]')
ylabel('Moment [Nm]')
%% Von mises

maxMy = max(abs(My));               %maxv�rde My
maxMz = max(abs(Mz));               %maxv�rde Mz
maxM  = sqrt(maxMy^2 + maxMy^2);    %Sammansatt b�jmoment enligt pythagoras

maxTy = max(abs(Ty));               %maxv�rde Ty
maxTz = max(abs(Tz));               %maxv�rde Tz
maxT  = sqrt(maxTy^2 + maxTz^2);    %Sammansatt tv�rkraft enligt pythagoras

Sigma_s = 310*10^6;            %Str�ckgr�ns 310Mpa f�r St�l SS165001 
Sigma_eff = Sigma_s/ns;        %Effektivsp�nning (ns given s�kerhetsfaktor)

syms d Wb Wv Sigma_max T_max            %DIMENSIONERING  

qn1 = (pi*d^3)/32 == Wb;            
qn2 = (pi*d^3)/16 == Wv;            
%qn3 = Sigma_max - maxM/Wb == 0;
%qn4 = T_max - Mv/Wv == 0;               
%qn5 = Sigma_eff - sqrt(Sigma_max^2+3*(T_max)^2) == 0; %Effektivsp�nning

s = solve([qn1, qn2], [d, Wb, Wv, Sigma_max, T_max]);
d = getfield(s,'d'); 
d = vpa(d); 
fprintf(' d = %g m\n ', d)

Mtot = sqrt(My.^2 + Mz.^2)';
figure(7)
plot(ix,Mtot)
title(['Totala Momentet Lastfall ' [LF]])
xlabel('L�ngd [m]')
ylabel('Moment [Nm]')

if xi < b1
       z = d/2;
   elseif xi < L-b1
       z = D/2;
   else
       z = d/2;
end

Wv = (pi*z^3)/16            %B�jmotst�nd  Se FS 6.9
Wb = (pi*z^3)/32            %Vridmotst�nd Se FS 6.78

Mtot = sqrt(Mxy.^2 + Mxz.^2) %vektoraddera och s�tt i en vektor
Sigma_max = Mtot/Wb
t_max = Mx/Wv %Vridmoment/Vridmotst�nd
Sigma_vm = []
Sigma_vm = sqrt(abs((Sigma_max).^2+3*(t_max)^2));
figure(8)
plot(Sigma_vm)
hold on
plot(Mtot)