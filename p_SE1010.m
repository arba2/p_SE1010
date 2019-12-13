%% Lastfall A - Konstant fart v
clear all, close all, clc, p_SE1010_MP; LF = 'A'; disp(['Lastfall ' LF])
a = 0;                     %Acceleration
v = v/2;                   %Fart
%% Lastfall B - Konstant accelaration a fr�n l�g fart
clear all, close all, clc, p_SE1010_MP; LF = 'B'; disp(['Lastfall ' LF])
a = a1;                    %Acceleration
v = v/2;                   %Fart
%% Lastfall C - Konstant retardation a fr�n maxfart
clear all, close all, clc, p_SE1010_MP; LF = 'C'; disp(['Lastfall ' LF])
a = a2;                    %Acceleration
v = v;                     %Fart
%% Lastfall D - Kurvtagning
clear all, close all, clc, p_SE1010_MP; LF = 'D'; disp(['Lastfall ' LF])
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
Vbi = FD/2;              %(FD/(Hbi + Hby))*Hbi   %OM LIVET VORE R�TTVIST
Vby = Vbi;               %(FD/(Hbi + Hby))*Hby

syms Hli Fli Fly Vli Vly
eq21 = Hbi - Hli + Hby == 0;                                                        %Kraftj�mvikt i x-led
eq22 = Vbi - Fli + Fk - Fly + Vby == 0;                                             %Kraftj�mvikt i y-led
eq23 = Vbi - Vli - Vb - Vly + Vby == 0;                                             %Kraftj�mvikt i z-led
%eq24 = - Vbi*rh - Mb*rb + Fk*rd - Vby*rh == 0;   �VERFL�DIG                        %Momentj�mvikt kring x
eq25 = Hbi*rh + Vli*b1 + Vb*(L/2 - bb) + Vly*(L - b1) + Hby*rh - Vby*L == 0;        %Momentj�mvikt kring y
eq26 = Vbi*0 - Fli*b1 + Fk*(L/2 + bd) - Fly*(L - b1) + Vby*L == 0;                  %Momentj�mvikt kring z-axeln

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
disp('Snittning')
%SP = [b1 L/2-bb L/2+bd L-b1 L];         %Snittpunkter

%0 <= x < b1
Tz1 = @() -Vbi;                          %Tv�rsp�nning Z-axeln
M1 = @(x) -Vbi * x - Hbi * rh;          %Moment XZ-planet

%b1 <= x < L/2-bb
Tz2 = @() Tz1() + Vli;                    %Tv�rsp�nning Z-axeln
My2 = @(x) M1(x) + Vli * (x-b1);         %Moment XZ-planet

%L/2-bb <= x < L-b1
Tz3 = @() Tz2() - Vb;                    %Tv�rsp�nning Z-axeln
My3 = @(x) My2(x) - Vb * (x-(L/2-bb));    %Moment XZ-planet

%L-b1 <= x < L
Tz4 = @() Tz3() + Vly;                    %Tv�rsp�nning Z-axeln
My4 = @(x) My3(x) + Vly * (x-(L-b1));     %Moment XZ-planet

%------------------------------------------------------
%0 <= x < b1
Ty1 = @() -FD/2;                         %Tv�rsp�nning Y-axeln
Mz1 = @(x) FD/2 * x;                    %Moment XY-planet

%b1 <= x < L/2+bd
Ty2 = @() Fli + Ty1();                    %Tv�rsp�nning Y-axeln
Mz2 = @(x) Mz1(x) - Fli * (x-b1);         %Moment XY-planet

%L/2-bd <= x < L-b1
Ty3 = @() -Fk + Ty2();                    %Tv�rsp�nning Y-axeln
Mz3 = @(x) Mz2(x) + Fk * (x-(L/2+bd));    %Moment XY-planet

%L-b1 <= x < L
Ty4 = @() Fly + Ty3();                    %Tv�rsp�nning Y-axeln
Mz4 = @(x) Mz3(x) - Fly * (x-(L-b1));     %Moment XY-planet

%------------------------------------------------------
%0 <= x < L/2-bb
Mx1 = @() -FD/2*rh;                      %Moment YZ-planet

%L/2-bb <= x < L/2+bd
Mx2 = @() -Vb*rb + Mx1();                %Moment YZ-planet

%L/2+bd <= x < L
Mx3 = @() Fk*rd + Mx2();                %Moment YZ-planet

%------------------------------------------------------
%0 <= x < b1
Nx1 = @() Hbi;                           %Normalkraft i X-axeln

%b1 <= x < L
Nx2 = @() -Hli + Nx1();                  %Normalkraft i X-axeln
%------------------------------------------------------
syms x

N(x)  = piecewise(0 <= x < b1, Nx1(), b1 <= x < L, Nx2());
Mx(x) = piecewise(0 <= x < L/2-bb, Mx1(), L/2-bb <= x < L/2+bd, Mx2(), L/2+bd <= x < L, Mx3());
Mz(x) = piecewise(0 <= x < b1, Mz1(x), b1 <= x < L/2+bd, Mz2(x), L/2-bd <= x < L-b1, Mz3(x), L-b1 <= x < L, Mz4(x));
Ty(x) = piecewise(0 <= x < b1, Ty1(), b1 <= x < L/2+bd, Ty2(), L/2-bd <= x < L-b1, Ty3(), L-b1 <= x < L, Ty4());
My(x) = piecewise(0 <= x < b1, M1(x), b1 <= x < L/2-bb, My2(x), L/2-bb <= x < L-b1, My3(x), L-b1 <= x < L, My4(x));
Tz(x) = piecewise(0 <= x < b1, Tz1(), b1 <= x < L/2-bb, Tz2(), L/2-bb <= x < L-b1, Tz3(), L-b1 <= x < L, Tz4());

clear z
z(x) = piecewise(0 < x <= b1, d/2, b1 < x < L-b1, D/2, L-b1 <= x < L, d/2);

%% ------------VON MISES------------
disp('von Mises')
%Wb = pi*z^3/32;                 %Se FS 6.9      %B�jmotst�nd
%Wv = pi*z^3/16;                 %Se FS 6.78     %Vridmotst�nd

I =  @(zt) pi*zt^4/4;            %Areatr�thetsmoment     Se FS 30.1.3
Wb = @(zt) I(zt)/abs(zt);        %B�jmotst�nd            Se G.L. S79
Wv = @(zt) pi*zt^3/2;            %Vridmotst�nd           Se FS 6.78

Aaxel = @(zt) zt^2*pi;           %Axel tv�rsnittsarea

Mtot = @(xt) sqrt(My(xt)^2 + Mz(xt)^2);                   %Sammansatt b�jmoment
Smax = @(xt) N(xt) / Aaxel(z(xt)) + Mtot(xt) / Wb(z(xt)); %Maxsp�nningen i axeln
Tmax = @(xt) Mx(xt) / Wv(z(xt));                          %Max skjuvsp�nning
VM   = @(xt) sqrt(Smax(xt)^2 + 3*Tmax(xt)^2);

%% Ber�kning av lokala sp�nningskoncentrationer vid �verg�ngar 
%...samt best�mning av axeldiameter D

SSmax = Smax./ns;    %Maximal till�ten sp�nning, inr�knat s�kerhetsfaktor ns.
D = 0.001;
d = 0.6*D;
KN = 1.45;                                          %Se F.S. 32.4 & 32.5
KM = 1.35;
KMx = 1.2;
SnomN  = @(x, d) double(4*N(x)/(pi*d^2));      %Nominella sp�nningar
SnomM  = @(x, d) double(32*Mtot(x)/(pi*d^3));
SnomMx = @(x, d) double(16*Mx(x)/(pi*d^3));

while true
    SmaxN1  = KN*SnomN(b1, d);
    SmaxM1  = KM*SnomM(b1, d);
    SmaxMx1 = KMx*SnomMx(b1, d);
    
    SmaxN2  = KN*SnomN(L-b1, d);
    SmaxM2  = KM*SnomM(L-b1, d);
    SmaxMx2 = KMx*SnomMx(L-b1, d);
    if SSmax > SmaxN1 && SSmax > SmaxM1 && SSmax > SmaxMx1 && SSmax > SmaxN2 && SSmax > SmaxM2 && SSmax > SmaxMx2
        break
    end
    D = D+0.001;
    d = 0.6*D;
end
disp([newline 'Slutgiltig n�dv�ndig diameter D f�r lastfall ' lastfall ': ' num2str(D*1000) ' mm.'])
disp(['d: ' num2str(d*1000) ' mm.'])
disp(['V�rde f�r grafer anges i Projekt_SE1010_variabler!' ])
fprintf('\nStr�ckgr�ns av material: %g Pa\n\nNormalsp�nning vid b1:   %g Pa\nB�jmoment vid b1:        %g Pa\nVridmoment vid b1:       %g Pa\nNormalsp�nning vid L-b1: %g Pa\nB�jmoment vid L-b1:      %g Pa\nVridmoment vid L-b1:     %g Pa\n\n ',SSmax,SmaxN1,SmaxM1,SmaxMx1,SmaxN2,SmaxM2,SmaxMx2)

%% UTMATTNING; reduktion av utmattningsdata figur 163 GH s.255 Ex. 43 s.259 GH
% Materialval Tab. 33.1 s.386 i FS , #7 SIS-141650-01
Rm  = 590; %MPa Brottgr�ns
Su  = 200; %+-200MPa sigma u, betecknar utmattningsgr�nsen vid v�xlande drag/tryck
Sup = 180; %180+-180MPa sigma up, betecknar utmattningsgr�nsen vid pulserande drag eller tryck
Subp = 240; %240+-240MPa sigma ubp, betecknar utmattningsgr�nsen vid pulserande b�jning
Ss  = 310; %sigma s >310 MPa Str�ckgr�ns
Ktd = 1.65;       %Sp�nningskoncentrationsfaktorn, drag figur 159b GH s.252
Ktb = 1.45;       %Sp�nningskoncentrationsfaktorn, b�j figur 159c GH s.252
q   = 0.85;        %K�lk�nslighetsfaktorn Figur 160 GH s.252

lambda = 1;%Teknologisk dimensionsfaktor. Axel ej gjuten figur 163 GH s.255
Kfd = 1 + q*(Ktb-1);  %Anvisningsfaktorn drag Ekv.(13-12) GH s.252
Kfb = 1 + q*(Ktb-1);  %Anvisningsfaktorn b�j Ekv.(13-12) GH s.252
Kd = 1;             %Geometriska volymsfaktorn figur 163 GH s.255
Kr = 1/0.975;   %Ra = 0.8my m enligt Ytfinhet wiki turning, figur 162b GH s.254

redfd = lambda/(Kfd*Kd*Kr); %reduktionsfaktor drag
redfb = lambda/(Kfb*Kd*Kr); %reduktionsfaktor drag Blir v�ldigt lika
redf = redfd
%Kdd = 1/0.96; %Geometriska volymsfaktorn figur 161 GH s.253
%KdD = 1/0.925

Haigh(x) = piecewise(0<=x<Sup, Su+((Sup-Su)/Sup)*x, Sup<=x<=Rm, Sup+(Sup/(Rm-Sup))*Sup+(Sup/(Sup-Rm))*x);
redSu = Su*redf;
redSup = Sup*redf;
redHaigh(x) = piecewise(0<=x<Sup, redSu+((redSup-redSu)/Sup)*x, Sup<=x<=Rm, redSup+(redSup/(Rm-Sup))*Sup+(redSup/(Sup-Rm))*x);
LSs = @(x) Ss-x %Linjen sigma s

fplot(Haigh,[0 Rm])
xlabel('Sigma a [MPa]')
ylabel('Sigma m [MPa]')
title('Haighdiagram')
grid on
axis equal
hold on
fplot(redHaigh,[0 Rm])

fplot(LSs, [0 Ss])
legend('Haighdiagram','Reducerat Haighdiagram (Arbetslinje)','Sigma s')

%Mittsp�nning och amplitud GH s.245 Nominella sp�nningar
Sm = 0; %(Smax + smin)/2 %mittsp�nning sigma m = 0 vid rent v�xlande sp�nning
Sa = SnomM; %S(Smax - Smin)/2 %sp�nningsamplitud sigma a
R = -1; %Sp�nningsf�rh�llande 
% OB = n*OA Belastningslinjen ligger p� y-axeln
%% PLOTTAR
close all

figure('Name',['Moment runt Y - Lastfall ' LF]);
fplot(My,[0 L]);
%legend("Moment XZ-planet runt Y")
xlabel('Position l�ngs X-axeln [m]')
ylabel('Moment [Nm]')
title(['Moment runt Y - Lastfall ' LF])
grid on

figure('Name',['Moment runt Z - Lastfall ' LF]);
fplot(Mz,[0 L]);
%legend("Moment XY-planet runt Z")
xlabel('Position l�ngs X-axeln [m]')
ylabel('Moment [Nm]')
title(['Moment  runt Z - Lastfall ' LF])
grid on

figure('Name',['Moment runt X - Lastfall ' LF]);
fplot(Mx,[0 L]);
%legend("Moment YZ-planet runt X")
xlabel('Position l�ngs X-axeln [m]')
ylabel('Moment [Nm]')
title(['Moment runt X - Lastfall ' LF])
grid on

figure('Name',['Tv�rsp�nning Z - Lastfall ' LF]);
fplot(Tz1,[0 L]);
%legend("Tv�rsp�nning Z-axeln")
xlabel('Position l�ngs X-axeln [m]')
ylabel('Tryck [Pa]')
title(['Tv�rsp�nning Z - Lastfall ' LF])
grid on

figure('Name',['Tv�rsp�nning Y - Lastfall ' LF]);
fplot(Ty,[0 L]);
%legend("Tv�rsp�nning Y-axeln")
xlabel('Position l�ngs X-axeln [m]')
ylabel('Tryck [Pa]')
title(['Tv�rsp�nning Y - Lastfall ' LF])
grid on

figure('Name',['Normalkraft X - Lastfall ' LF]);
fplot(N,[0 L]);
%legend("Normalkraft mot YZ-planet")
xlabel('Position l�ngs X-axeln [m]')
ylabel('Kraft [N]')
title(['Normalkraft X - Lastfall ' LF])
grid on

VMt = [];
for i = 0:0.01:L
   VMt = [VMt VM(i)]; 
end
figure('Name',['Von Mises - Lastfall ' LF]);
plot(0:0.01:L,VMt)
%fplot(VM,[0 L]);
%legend("Von Mises")
xlabel('Position l�ngs X-axeln [m]')
ylabel('Moment [Nm] / Kraft [N]')
title(['Von Mises - Lastfall ' LF])
grid on