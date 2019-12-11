% Givna variabler
    A   = 0.6;              %m2         Fortonets frontarea
    a1  = 6;                %m/s2       Max acceleration
    a2  = -15;              %m/s2       Max retardation    
    b1  = 0.15;             %m          Hjullagerposition
    bb  = 0.3;              %m          Bromsskivans position
    bd  = 0.1;              %m          Drevets position
    c   = 0.3;  %0.4;              %           Luftmotståndskoefficient
    db  = 0.3;  %0.4;              %m          Avstånd från bakaxel till tyngdpunkt
    df  = 0.8;  %1;                %m          Avstånd från framaxel till tyngpunkt
    dh  = 0.3;              %m          Hjuldiameter
    h   = 0.3;  %0.6;              %m          Tyngdpunktens vertikala position
    h1  = 0.3;  %0.2;              %m          Avstånd mellan tyngdpunkt och luftmotståndets verkningslinje
    L   = 1.1;    %1.15;             %m          Bakaxellängd
    m   = 120; %200;              %kg         Fordonets totalvikt
    ns  = 3;                %           Säkerhetsfaktor mot plastisk deformation
    nu  = 1.8;              %           Säkerhetsfaktor mot utmattning
    R   = 20;               %m          Kurvradie (vid 0,5v)
    rh  = dh/2;             %m          Hjulradie
    rb  = 0.03;             %m          Bromsbackarnas position
    rd  = 0.03;             %m          Drevets radie
    v   = 100/3.6; %150/3.6;          %m/s        Max fart (rakt fram) 
    p1  = 0.005;            %m          Kälradie
    p2  = p1;               %m          Kälradie
    % D söks
    t   = 1;                %s          Förfluten tid sedan påbörjad acceleration
% Övriga variabler
    g   = 9.82;             %m/s2       Tyngdacceleration
    pluft = 1.21;  %1.2041;         %kg/m3      Densitet luft
    D  = 0.044;
    Hc = 0;                  %Centrifugalkraften
%{
            Mer
            FL = Luftmotstånd
            FD = Friktionskraft mellan bägge makhjul och väg

            INDEXERING 2.0
            V = Vertikal    (upp-ner) ((tidigare R))
            H = Horisontell (höger-vänster)
            F =             (fram-bak)
            f = fram
            b = bak
            h = hjul
            i = inre
            y = yttre
    
            KOORDINATER
            Origo = axelns centrumpunkt
            x = ut till höger
            y = fram
            z = upp
            
    
            *********************GAMMALLT******************************
            KRAFTER SOM ANGRIPER I AXELN
            RH = Hjulets bärande kraft (trycker upp/ner)
            FH = Hjulets
            MH = Hjulets vridande moment på axeln
            RLager = Lagrens tyngande kraft på axeln
            FLager = Lagrens förskjutande kraft i y-led (trycker fram/bak)
            Fk = Kedjekraft
            Md = Momentet från drevet
            Fb = Bromskraft
            Md = Momentet av bromsskivan
            Hf = Horisontell kraft på framhjul
            Hbi = Horisontell kraft inre bakhjul
            Hby = Horisontell kraft yttre bakhjul
    
            Rf = reaktionskraft  framhjul
            Rb = reaktionskraft bägge bakhjul
    %}