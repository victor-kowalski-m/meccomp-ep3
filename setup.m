clear
clc
close all

dx = 3; % discretização escolhida

%%% Pontos do desenho
pontos.A.coor = [0  0 ]; pontos.A.con = 1; % A
pontos.B.coor = [5  0 ]; pontos.B.con = 0; % B
pontos.C.coor = [9  0 ]; pontos.C.con = 0; % C
pontos.D.coor = [18 0 ]; pontos.D.con = 0; % D
pontos.E.coor = [36 0 ]; pontos.E.con = 0; % E
pontos.F.coor = [36 18]; pontos.F.con = 0; % F
pontos.G.coor = [36 -5]; pontos.G.con = 0; % G
pontos.H.coor = [65  0]; pontos.H.con = 0; % H

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
    pontos.(ps(i+1)).con = nos;
end
vigaHor = 1:size(con, 1);

%%% Insere viga vertical nas matrizes de coord. e conec.
for p=['G', 'F']
    [pIn, NpIn] = calcula_pontos_internos( ...
        pontos.E.coor, pontos.(p).coor, dx);
    pIn = pIn(2:end, :); NpIn = NpIn - 1;
    coor = [coor; pIn];
    con = [con; [pontos.E.con nos+1]]; 
    con = [con; nos+[(1:NpIn-1)', (2:NpIn)']];
    nos = nos + NpIn;
    pontos.(p).con = nos;
end
vigaVer = vigaHor(2)+1:size(con, 1);

%%% Insere treliças nas matrizes de coord. e conec.
for p = ['B', 'D', 'H']
    con(end+1, :) = [pontos.(p).con pontos.F.con];
end
trelicas = vigaVer(2)+1:size(con, 1);

% plot_struct(coor, con)

%%% Quantidades de nós, elementos e G.L.
Nod = size(coor, 1);
Nel = size(con, 1);
Ngdl = 3*Nel;

%%% Props. dos elementos
data.E = 210e9*ones(Nel, 1);
data.rho = 7600*ones(Nel, 1);

data.A = zeros(Nel, 1);
data.A(vigaHor) = 1.8*0.9; 
data.A(vigaVer) = 0.9*0.9;
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
end


%%% Funções para cálculo das matrizes de rigidez, massa e transformação
Ke = @(E, A, I, L)...
    [ ...
     1  0  0 -1  0  0;
     0  0  0  0  0  0;
     0  0  0  0  0  0;
    -1  0  0  1  0  0;
     0  0  0  0  0  0;
    ] * E*A/L + ...
    [ ...
     0      0      0      0     0     0     ;
     0      12     6*L    0    -12    6*L   ;
     0      6*L    4*L^2  0    -6*L   2*L^2 ;
     0      0      0      0     0     0     ;
     0     -12    -6*L    0     12   -6*L   ;
     0      6*L    2*L^2  0     -6*L  4*L^2 ;
    ] * E*I/L^3 ;

Me = @(rho, A, L)...
    [ ...
     2  0  0  1  0  0;
     0  0  0  0  0  0;
     0  0  0  0  0  0;
     1  0  0  2  0  0;
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

T = @(theta)...
    [
     cos(theta) -sin(theta)  0          0           0           0;
     sin(theta)  cos(theta)  0          0           0           0;
     0           0           1          0           0           0;
     0           0           0          cos(theta) -sin(theta)  0;
     0           0           0          sin(theta)  cos(theta)  0;
     0           0           0          0           0           1
    ];

Kg = zeros(Ngdl);
Mg = zeros(Ngdl);




