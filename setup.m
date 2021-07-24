clear
clc
close all

dx = 1; % discretização escolhida

%%% Pontos do desenho
pontos.A.coor = [0  0 ]; pontos.A.nod = 1; % A
pontos.B.coor = [5  0 ]; pontos.B.nod = 0; % B
pontos.C.coor = [9  0 ]; pontos.C.nod = 0; % C
pontos.D.coor = [18 0 ]; pontos.D.nod = 0; % D
pontos.E.coor = [36 0 ]; pontos.E.nod = 0; % E
pontos.F.coor = [36 18]; pontos.F.nod = 0; % F
pontos.G.coor = [36 -5]; pontos.G.nod = 0; % G
pontos.H.coor = [65  0]; pontos.H.nod = 0; % H

coor = pontos.A.coor; con = []; nos = 1; % inicialização de matrizes de
    % coordenadas e conectividade

%%% Insere viga horizontal nas matrizes de coord. e conec.
ps = ['A', 'B', 'C', 'D', 'E', 'H'];
for i=1:length(ps)-1
    [pIn, NpIn] = calcula_pontos_internos( ...
        pontos.(ps(i)).coor, pontos.(ps(i+1)).coor, dx);
    pIn = pIn(2:end, :); NpIn = NpIn - 1;
    coor = [coor; pIn]; 
    con = [con; nos+[(0:NpIn-1)', (1:NpIn)']];
    nos = nos + NpIn;
    pontos.(ps(i+1)).nod = nos;
end
vigaHor = 1:size(con, 1);

%%% Insere viga vertical nas matrizes de coord. e conec.
for p=['G', 'F']
    [pIn, NpIn] = calcula_pontos_internos( ...
        pontos.E.coor, pontos.(p).coor, dx);
    pIn = pIn(2:end, :); NpIn = NpIn - 1;
    coor = [coor; pIn];
    con = [con; [pontos.E.nod nos+1]]; 
    con = [con; nos+[(1:NpIn-1)', (2:NpIn)']];
    nos = nos + NpIn;
    pontos.(p).nod = nos;
end
vigaVer = (vigaHor(end)+1):size(con, 1);

%%% Insere treliças nas matrizes de coord. e conec.
for p = ['B', 'D', 'H']
    con(end+1, :) = [pontos.(p).nod pontos.F.nod];
end
trelicas = (vigaVer(end)+1):size(con, 1);

plot_struct(coor, con)

%%% Quantidades de nós, elementos e G.L.
Nod = size(coor, 1);
Nel = size(con, 1);
Ngdl = 3*Nod;

%%% Props. dos elementos
data.E = 210e9*ones(Nel, 1);
data.rho = 7600*ones(Nel, 1);

data.A = zeros(Nel, 1);
data.A(vigaHor) = 0.9*0.9; 
data.A(vigaVer) = 1.8*0.9;
data.A(trelicas) = pi*0.05^2/4;

data.L = zeros(Nel, 1);
data.Q = zeros(Nel, 1);

for e = 1:Nel
    
    x1 = coor(con(e, 1), 1);
    x2 = coor(con(e, 2), 1);
    
    y1 = coor(con(e, 1), 2);
    y2 = coor(con(e, 2), 2);
    
    data.L(e) = sqrt((x2-x1)^2 + (y2-y1)^2);
    data.Q(e) = atan2(y2-y1, x2-x1)*180/pi;
    data.I(e) = data.rho(e)*data.A(e)*data.L(e)^3/12;
    
end


%%% Funções para cálculo das matrizes de rigidez, massa e transformação
Ke = @(E, A, I, L)...
    [ ...
     1  0  0 -1  0  0;
     0  0  0  0  0  0;
     0  0  0  0  0  0;
    -1  0  0  1  0  0;
     0  0  0  0  0  0;
     0  0  0  0  0  0;
    ] * E*A/L + ...
    [ ...
     0      0      0      0     0     0     ;
     0      12     6*L    0    -12    6*L   ;
     0      6*L    4*L^2  0    -6*L   2*L^2 ;
     0      0      0      0     0     0     ;
     0     -12    -6*L    0     12   -6*L   ;
     0      6*L    2*L^2  0    -6*L   4*L^2 ;
    ] * E*I/(L^3) ;

Me = @(rho, A, L)...
    [ ...
     2  0  0  1  0  0;
     0  0  0  0  0  0;
     0  0  0  0  0  0;
     1  0  0  2  0  0;
     0  0  0  0  0  0;
     0  0  0  0  0  0;
    ] * rho*A*L/6 + ...
    [ ...
     0       0       0       0      0      0     ;
     0       156     22*L    0      54    -13*L   ;
     0       22*L    4*L^2   0      13*L  -3*L^2 ;
     0       0       0       0      0      0     ;
     0       54      13*L    0      156   -22*L   ;
     0      -13*L   -3*L^2   0     -22*L   4*L^2 ;
    ] * rho*A*L/420 ;

Te = @(theta)...
    [
     cosd(theta)  sind(theta)  0          0           0            0;
    -sind(theta)  cosd(theta)  0          0           0            0;
     0            0            1          0           0            0;
     0            0            0          cosd(theta) sind(theta)  0;
     0            0            0         -sind(theta) cosd(theta)  0;
     0            0            0          0            0           1
    ];

Kg = zeros(Ngdl);
Mg = zeros(Ngdl);

for e = 1:Nel
    
    ke = Ke(data.E(e), data.A(e), data.I(e), data.L(e));
    me = Me(data.rho(e), data.A(e), data.L(e));
    
    T = Te(data.Q(e));
    ke = T'*ke*T;
    me = T'*me*T;

    nod1 = con(e, 1);
    nod2 = con(e, 2);
    
    Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = ...
        Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + ke(1:3, 1:3);
    Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = ...
        Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + ke(1:3, 4:6);
    Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = ...
        Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + ke(4:6, 1:3); 
    Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = ...
        Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + ke(4:6, 4:6);
    
    Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = ...
        Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + me(1:3, 1:3);
    Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = ...
        Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + me(1:3, 4:6);
    Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = ...
        Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + me(4:6, 1:3); 
    Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = ...
        Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + me(4:6, 4:6);
        
end

Kgm = Kg;
Mgm = Mg;

% uE = 3*pontos.E.nod-2;
% wE = 3*pontos.E.nod-1;
% Kgm(:, uE) = 0; Kgm(uE, :) = 0; Kgm(uE, uE) = 1;
% Kgm(:, wE) = 0; Kgm(wE, :) = 0; Kgm(wE, wE) = 1;
% Mgm(:, uE) = 0; Mgm(uE, :) = 0; Mgm(uE, uE) = 1;
% Mgm(:, wE) = 0; Mgm(wE, :) = 0; Mgm(wE, wE) = 1; 

uH = 3*pontos.H.nod-2;
wH = 3*pontos.H.nod-1;
Kgm(:, uH) = 0; Kgm(uH, :) = 0; Kgm(uH, uH) = 1;
Kgm(:, wH) = 0; Kgm(wH, :) = 0; Kgm(wH, wH) = 1;
Mgm(:, uH) = 0; Mgm(uH, :) = 0; Mgm(uH, uH) = 1;
Mgm(:, wH) = 0; Mgm(wH, :) = 0; Mgm(wH, wH) = 1;

uG = 3*pontos.G.nod-2;
wG = 3*pontos.G.nod-1;
phiG = 3*pontos.G.nod;
Kgm(:, uG) = 0; Kgm(uG, :) = 0; Kgm(uG, uG) = 1;
Kgm(:, wG) = 0; Kgm(wG, :) = 0; Kgm(wG, wG) = 1;
Kgm(:, phiG) = 0; Kgm(phiG, :) = 0; Kgm(phiG, phiG) = 1;
Mgm(:, uG) = 0; Mgm(uG, :) = 0; Mgm(uG, uG) = 1;
Mgm(:, wG) = 0; Mgm(wG, :) = 0; Mgm(wG, wG) = 1;
Mgm(:, phiG) = 0; Mgm(phiG, :) = 0; Mgm(phiG, phiG) = 1; 

A = Mgm\Kgm;
[vec, val] = eig(A);
val = sqrt(diag(val))/(2*pi);

