% Givna variabler
    A   = 0.6;              %m2         Fortonets frontarea
    a1  = 6;                %m/s2       Max acceleration
    a2  = -15;              %m/s2       Max retardation    
    b1  = 0.15;             %m          Hjullagerposition
    bb  = 0.3;              %m          Bromsskivans position
    bd  = 0.1;              %m          Drevets position
    c   = 0.3;  %0.4;              %           Luftmotst�ndskoefficient
    db  = 0.3;  %0.4;              %m          Avst�nd fr�n bakaxel till tyngdpunkt
    df  = 0.8;  %1;                %m          Avst�nd fr�n framaxel till tyngpunkt
    dh  = 0.3;              %m          Hjuldiameter
    h   = 0.3;  %0.6;              %m          Tyngdpunktens vertikala position
    h1  = 0.3;  %0.2;              %m          Avst�nd mellan tyngdpunkt och luftmotst�ndets verkningslinje
    L   = 1.1;    %1.15;             %m          Bakaxell�ngd
    m   = 120; %200;              %kg         Fordonets totalvikt
    ns  = 3;                %           S�kerhetsfaktor mot plastisk deformation
    nu  = 1.8;              %           S�kerhetsfaktor mot utmattning
    R   = 20;               %m          Kurvradie (vid 0,5v)
    rh  = dh/2;             %m          Hjulradie
    rb  = 0.03;             %m          Bromsbackarnas position
    rd  = 0.03;             %m          Drevets radie
    v   = 100/3.6; %150/3.6;          %m/s        Max fart (rakt fram) 
    p1  = 0.005;            %m          K�lradie
    p2  = p1;               %m          K�lradie
    % D s�ks
    t   = 1;                %s          F�rfluten tid sedan p�b�rjad acceleration
% �vriga variabler
    g   = 9.82;             %m/s2       Tyngdacceleration
    pluft = 1.21;  %1.2041;         %kg/m3      Densitet luft
    D  = 0.044;
    Hc = 0;                  %Centrifugalkraften
%{
            Mer
            FL = Luftmotst�nd
            FD = Friktionskraft mellan b�gge makhjul och v�g

            INDEXERING 2.0
            V = Vertikal    (upp-ner) ((tidigare R))
            H = Horisontell (h�ger-v�nster)
            F =             (fram-bak)
            f = fram
            b = bak
            h = hjul
            i = inre
            y = yttre
    
            KOORDINATER
            Origo = axelns centrumpunkt
            x = ut till h�ger
            y = fram
            z = upp
            
    
            *********************GAMMALLT******************************
            KRAFTER SOM ANGRIPER I AXELN
            RH = Hjulets b�rande kraft (trycker upp/ner)
            FH = Hjulets
            MH = Hjulets vridande moment p� axeln
            RLager = Lagrens tyngande kraft p� axeln
            FLager = Lagrens f�rskjutande kraft i y-led (trycker fram/bak)
            Fk = Kedjekraft
            Md = Momentet fr�n drevet
            Fb = Bromskraft
            Md = Momentet av bromsskivan
            Hf = Horisontell kraft p� framhjul
            Hbi = Horisontell kraft inre bakhjul
            Hby = Horisontell kraft yttre bakhjul
    
            Rf = reaktionskraft  framhjul
            Rb = reaktionskraft b�gge bakhjul
    %}