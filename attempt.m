
ix = [0, b1, L/2 - bb, L/2 + bd, L - b1, L];   %Snittpunkter

x<ix(2)
fs = {Nx1, Nx2, Nx3, Nx4, Nx5                  %FunktionsSparare
      Ty1, Ty2...} 
syms x
Tz = piecewise(sp(1)<=x<sp(2), Vbi              ,sp(2)<=x<sp(3),                 

My = piecewise(sp(1)<=x<sp(2), -Hbi*rh - Vbi*x  ,sp(2)<=x<sp(3), 
Ty = piecewise(sp(1)<=x<sp(2), Fbi              ,sp(2)<=x<sp(3), 
Mz = piecewise(sp(1)<=x<sp(2), Fbi*x            ,sp(2)<=x<sp(3), 
Nx = piecewise(sp(1)<=x<sp(2), - Hbi + Hli - Hby,sp(2)<=x<sp(3), 
Mx = piecewise(sp(1)<=x<sp(2), FD*rh            ,sp(2)<=x<sp(3), 
fplot(y)
%xz T&M
e1 = 1 %elementet av sp och intervallets nummer
e2 = 2
x1 = sp(e1)
x2 = sp(e2)
while e2<5
for x = ix(1):h:ix(2)
    Nx = @(x) - Hbi + Hli - Hby
    Ty = @(x) Fbi
    Tz = @(x) Vbi
    Mx = @(x) FD*rh
    My = @(x) -Hbi*rh - Vbi*x
    Mz = @(x) Fbi*x
    
for x = ix(2):h:ix(3)    
    Ty = @(x) Fbi - Fli
    Tz = @(x) Vbi - Vli
    My = @(x) -Hbi*rh - Vbi*x + Vli*(x - b1)
    Mz = @(x) Fbi*x + Fli*(x - b1)
end
    
    
    e1 = e1 + 1
    e2 = e2 + 1
    
 if Snitt(1,i) < b1
        z = d/2;
    elseif Snitt(1,i) < L-b1
        z = D/2;
    else
        z = d/2;
    end
    Wb = pi*z^3/32;                 %Se FS 6.9
    Wv = pi*z^3/16;                 %Se FS 6.78
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

