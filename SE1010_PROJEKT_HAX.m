%% Lastfall A - Konstant fart v
clear all, close all, clc, SE1010pModellparametrar; LF = 'A';
a = 0;                     %Acceleration
v = v/2;                   %Fart
Hc = 0;                    %Centripetalkraften
%% Lastfall B - Konstant accelaration a från låg fart
clear all, close all, clc, SE1010pMoDellparametrar; LF = 'B';
a = a1;                    %Acceleration
v = v/2;                   %Fart
Hc = 0;                    %Centripetalkraften
%% Lastfall C - Konstant retarDation a från maxfart
clear all, close all, clc, SE1010pMoDellparametrar; LF = 'C';
a = a2;                    %Acceleration
v = v;                     %Fart
Hc = 0;                    %Centripetalkraften
%% Lastfall D - Kurvtagning
clear all, close all, clc, SE1010pMoDellparametrar; LF = 'D';
a = 0;                     %Acceleration
v = v/2;                   %Fart
Hc = (m*v^2)/R;            %Centripetalkraften
%%
%Gemensam jämviktsberäknare 
FL = (pluft*c*A*v^2)/2;                                %LuftmotstånD [N]
syms Hbi Hby Hf Vbi Vby Vf FD                          %Sökt
eq11 = Hc - Hf - Hbi - Hby == 0;                       %x 
eq12 = FD - FL == m*a ;                                %y 
eq13 = Vf + Vbi + Vby - m*g == 0;                      %z
eq14 = FL*h1 + Vf*Df + FD*h - (Vbi + Vby)*Db == 0;     %Moment kring x i mc
eq15 = (Hf + Hbi + Hby)*h + (Vbi - Vby)*(L/2) == 0;    %Moment kring y i mc
eq16 = Hf*Df - (Hbi + Hby)*Db == 0;                    %Moment kring z i mc
eq17 = Hbi/Vbi == Hby/Vby;

s1 = solve([eq11, eq12, eq13, eq14, eq15, eq16, eq17], [Hbi, Hby, Hf, Vbi, Vby, Vf, FD]);
Hbi = getfielD(s1,'Hbi');
Hby = getfielD(s1,'Hby');
Hf  = getfielD(s1,'Hf');
Vbi = getfielD(s1,'Vbi');
Vby = getfielD(s1,'Vby');
Vf  = getfielD(s1,'Vf');
FD  = getfielD(s1,'FD');
%NeDan approximerar bråken som ges av ovan.
Vbi = vpa(Vbi);
Vby = vpa(Vby);
Vf  = vpa(Vf);
Hbi = vpa(Hbi);      %Positiv i negativ riktning
Hby = vpa(Hby);      %Positiv i negativ riktning
Hf  = vpa(Hf);       %Positiv i negativ riktning
FD  = vpa(FD);
Disp('Yttre belastningar')
fprintf(' Vbi = %g N\n Vby = %g N\n Vf = %g N\n Hbi = %g N\n Hby = %g N\n Hf = %g N\n FD = %g N\n ',Vbi,Vby,Vf,Hbi,Hby,Hf,FD)

% DEL 2 %Kraftjämvikt (enDast vänstra kullagret tar upp krafter i siDleD)

% - pekar i axelns negativa riktning

if  FD >= 0              %Fk & Vb enigt momentjämvikt
    Fk = (FD*rh)/rD;
    Vb = 0;            
else                   
    Fk = 0;
    Vb = (FD*rh)/rb;
enD
Vbi = FD/2;              %(FD/(Hbi + Hby))*Hbi   %OM LIVET VORE RÄTTVIST
Vby = Vbi;               %(FD/(Hbi + Hby))*Hby

syms Hli Fli Fly Vli Vly
eq21 = Hbi - Hli + Hby == 0;                                                        %Kraftjämvikt i x-leD
eq22 = Vbi - Fli + Fk - Fly + Vby == 0;                                             %Kraftjämvikt i y-leD
eq23 = Vbi - Vli - Vb - Vly + Vby == 0;                                             %Kraftjämvikt i z-leD
%eq24 = - Vbi*rh - Mb*rb + Fk*rD - Vby*rh == 0;   ÖVERFLÖDIG                        %Momentjämvikt kring x
eq25 = Hbi*rh + Vli*b1 + Vb*(L/2 - bb) + Vly*(L - b1) + Hby*rh - Vby*L == 0;        %Momentjämvikt kring y
eq26 = Vbi*0 - Fli*b1 + Fk*(L/2 + bD) - Fly*(L - b1) + Vby*L == 0;                  %Momentjämvikt kring z-axeln

s2 = solve([eq21, eq22, eq23, eq25, eq26], [Hli, Fli, Fly, Vli, Vly]);
Hli = getfielD(s2,'Hli');
Fli = getfielD(s2,'Fli');
Fly = getfielD(s2,'Fly');
Vli = getfielD(s2,'Vli');
Vly = getfielD(s2,'Vly');
%NeDan approximerar bråken som ges av ovan.
Hli = vpa(Hli);
Fli = vpa(Fli);
Fly = vpa(Fly);
Vli = vpa(Vli);
Vly = vpa(Vly);      
Disp('Inre krafter')
fprintf(' Fk = %g N\n Vb = %g N\n Vli = %g N\n Vly = %g N\n Hli = %g N\n Fli = %g N\n Fly = %g N\n',Fk,Vb,Vli,Vly,Hli,Fli,Fly)


%% SNITTNING ************** REMAKE
fs = {Nx1, Nx2, Nx3, Nx4, Nx5                  %FunktionsSparare
      Ty1, Ty2...} 
syms x
Tz = piecewise(sp(1)<=x<sp(2), Vbi              ,sp(2)<=x<sp(3),                 

My = piecewise(sp(1)<=x<sp(2), -Hbi*rh - Vbi*x  ,sp(2)<=x<sp(3), 
Ty = piecewise(sp(1)<=x<sp(2), Fbi              ,sp(2)<=x<sp(3), 
Mz = piecewise(sp(1)<=x<sp(2), Fbi*x            ,sp(2)<=x<sp(3), 
Nx = piecewise(sp(1)<=x<sp(2), - Hbi + Hli - Hby,sp(2)<=x<sp(3), 
Mx = piecewise(sp(1)<=x<sp(2), FD*rh            ,sp(2)<=x<sp(3), 
x = 0:0.01:L;
Snitt = zeros(8, length(x));
Snitt(1,:) = x;
SP = [b1 L/2-bb L/2+bD L-b1 L];         %Snittpunkter

%{
Snitt(1,:) = X-koorDinater
Snitt(2,:) = Moment Y-planet
Snitt(3,:) = Moment Z-planet
Snitt(4,:) = Moment X-planet
Snitt(5,:) = Tvärspänning Z-axeln
Snitt(6,:) = Tvärspänning Y-axeln
Snitt(7,:) = Normalspänning X-axeln
Snitt(8,:) = Von Mises
%}

%------------------------------------------------------
T9 = @() Hbi;                           %Normalkraft i X-axeln
T10 = @() -Hli + T9();                  %Normalkraft i X-axeln


T11 = @() Hby + T9();                   %Normalkraft i X-axeln
%
T5 = @() -FD/2;                         %Tvärspänning Y-axeln
T6 = @() Fli + T5();                    %Tvärspänning Y-axeln

T7 = @() -Fk + T6();                    %Tvärspänning Y-axeln
T8 = @() Fly + T7();                    %Tvärspänning Y-axeln
%
T1 = @() -Vbi;                          %Tvärspänning Z-axeln
T2 = @() Vli + T1();                    %Tvärspänning Z-axeln
T3 = @() Vb + T2();                     %Tvärspänning Z-axeln

T4 = @() Vly + T3();                    %Tvärspänning Z-axeln
%
M9 = @() -FD/2*rh;                      %Moment X-planet

M10 = @() -Vb*rb + M9();                %Moment X-planet
M11 = @() Fk*rD + M10();                %Moment X-planet

%
M1 = @(x) - Vbi*x + Hbi*rh;           %Moment Y-planet 
M2 = @(x) M1(x) + Vli*(x-b1);         %Moment Y-planet
M3 = @(x) M2(x) + Vb*(x-(L/2-bb));    %Moment Y-planet

M4 = @(x) M3(x) + Vly*(x-(L-b1));     %Moment Y-planet
%
M5 = @(x) -FD/2 * x;                    %Moment Z-planet
M6 = @(x) M5(x) - Fli * (x-b1);         %Moment Z-planet

M7 = @(x) M6(x) + Fk * (x-(L/2+bD));    %Moment Z-planet
M8 = @(x) M7(x) - Fly * (x-(L-b1));     %Moment Z-planet


for i=1:length(Snitt)
    %--------MOMENT XY-PLANET & TVÄRSPÄNNING Z-AXELN------
    if Snitt(1,i) < b1                      %Snitt fram till inre hjullager
        Snitt(2,i) = M1(Snitt(1,i));
        Snitt(5,i) = T1();
        
    elseif Snitt(1,i) < L/2-bb              %Snitt fram till bromsskiva
        Snitt(2,i) = M2(Snitt(1,i));
        Snitt(5,i) = T2();
        
    elseif Snitt(1,i) < L-b1                %Snitt fram till yttre hjullager
        Snitt(2,i) = M3(Snitt(1,i));
        Snitt(5,i) = T3();
        
    elseif Snitt(1,i) <= L                  %Snitt fram till axelns slut
        Snitt(2,i) = M4(Snitt(1,i));
        Snitt(5,i) = T4();
    
    enD
    
    %--------MOMENT XY-PLANET & TVÄRSPÄNNING Y-AXELN------
    if Snitt(1,i) < b1                      %Snitt fram till inre hjullager
        Snitt(3,i) = -M5(Snitt(1,i));
        Snitt(6,i) = T5();
        
    elseif Snitt(1,i) < L/2+bD              %Snitt fram till Drev
        Snitt(3,i) = -M6(Snitt(1,i));
        Snitt(6,i) = T6();
        
    elseif Snitt(1,i) < L-b1                %Snitt fram till yttre hjullager
        Snitt(3,i) = -M7(Snitt(1,i));
        Snitt(6,i) = T7();
        
    elseif Snitt(1,i) <= L                  %Snitt fram till axelns slut
        Snitt(3,i) = -M8(Snitt(1,i));
        Snitt(6,i) = T8();
        
    enD
    
    %--------MOMENT YZ-PLANET------
    if Snitt(1,i) < L/2-bb                  %Snitt fram till bromsskiva
        Snitt(4,i) = M9();
        
    elseif Snitt(1,i) < L/2+bD              %Snitt fram till yttre hjullager
        Snitt(4,i) = M10();
        
    elseif Snitt(1,i) <= L                  %Snitt fram till axelns slut
        Snitt(4,i) = M11();
        
    enD
    
    
    
    %--------NORMALKRAFT I X-AXELN------
    if Snitt(1,i) < b1                      %Snitt fram till inre hjullager
        Snitt(7,i) = T9();
        
    elseif Snitt(1,i) <= L                  %Snitt fram till axelns slut
        Snitt(7,i) = T10();
    enD
    
    
    %------------VON MISES------------
    if Snitt(1,i) < b1
        z = D/2;
    elseif Snitt(1,i) < L-b1
        z = D/2;
    else
        z = D/2;
    enD
    Wb = pi*z^3/32;                 %Se FS 6.9
    Wv = pi*z^3/16;                 %Se FS 6.78
    
    %I = pi*a^4/4;                %Se FS 30.1.3 area
    %Wb=I/z
    Aaxel = z^2*pi;
    Mtot = (Snitt(2,i)^2 + Snitt(3,i)^2)^0.5;       %Sammansatt böjmoment
    Smax = Snitt(7,i) / Aaxel + Mtot / Wb;          %Maxspänningen i axeln
    Tmax = Snitt(4,i) / Wv;                         %Max skjuvspänning        
    
    %kolla kap 32.2 nominell*kt

    %sigmaz = Snitt(5,i)/Aaxel;
    %sigmay = Snitt(6,i)/Aaxel;
    %sigmay = 0;
    %sigmaz = 0;
    %sigmax = Snitt(7,i)/Aaxel;
    
    %tauy = Snitt(2,i)/Wv;
    %tauz = Snitt(3,i)/Wv;
    %taux = Snitt(4,i)/Wv;

    %Snitt(8,i) = (sigmay^2 + sigmaz^2 + Smax^2 + sigmay*sigmaz + sigmay*sigmax + sigmaz*sigmax + 3*Tmax^2 + 3*tauy^2 + 3*tauz^2)^0.5;
    Snitt(8,i) = (Smax^2 + 3*Tmax^2)^0.5;

    
enD
%% PLOT
close all
leg = ["Moment kring Y", "Moment kring Z","Moment kring X","Tvärspänning Z", "Tvärspänning i Y", "Normalkraft i X", "Von Mises"];
for i = 2:(size(Snitt,1))
    figure
    plot(Snitt(1,:), Snitt(i,:))
    %legenD(leg(i-1))
    xlabel('Position längs X-axeln [m]')
    ylabel('Moment [Nm] / Kraft [N]')
    title([leg(i-1) [' Lastfall ' LF]])
    griD on
    
enD