function [coor, con, pontos, data, Kgm, Mgm, Fg, id_free, Ngdl]...
    = setup(dx)

% Pontos do desenho
pontos.A.coor = [0  0 ]; pontos.A.nod = 1; % A
pontos.B.coor = [5  0 ]; pontos.B.nod = 0; % B
pontos.C.coor = [9  0 ]; pontos.C.nod = 0; % C
pontos.D.coor = [18 0 ]; pontos.D.nod = 0; % D
pontos.E.coor = [36 0 ]; pontos.E.nod = 0; % E
pontos.F.coor = [36 18]; pontos.F.nod = 0; % F
pontos.G.coor = [36 -5]; pontos.G.nod = 0; % G
pontos.H.coor = [45  0]; pontos.H.nod = 0; % H
pontos.I.coor = [65  0]; pontos.I.nod = 0; % I

% Inicialização de matrizes de coordenadas e conectividade
coor = pontos.A.coor;
con = [];
nos = 1;

% Insere viga horizontal nas matrizes de coord. e conec.
if dx == 1
    ps = ['A', 'B', 'C', 'D', 'E', 'H','I'];
else
    ps = ['A', 'B', 'D', 'E', 'I'];
end
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

% Insere viga vertical nas matrizes de coord. e conec.
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
for p = ['B', 'D', 'I']
    con(end+1, :) = [pontos.(p).nod pontos.F.nod];
end
trelicas = (vigaVer(end)+1):size(con, 1);

% Quantidades de nós, elementos e G.L.
Nod = size(coor, 1);
Nel = size(con, 1);
Ngdl = 3*Nod;

% Modulo de Young e densidade dos elementos
data.E = 210e9*ones(Nel, 1);
data.rho = 7600*ones(Nel, 1);

% Areas dos elementos 
data.A = zeros(Nel, 1);
data.A(vigaHor) = 0.9*0.9;
data.A(vigaVer) = 1.8*0.9;
data.A(trelicas) = pi*0.05^2/4;

% Momentos de inércia dos elementos
data.I = zeros(Nel, 1);
data.I(vigaHor) = 0.9*0.9^3/12;
data.I(vigaVer) = 1.8^3*0.9/12;
data.I(trelicas) = pi*0.05^4/32;

% Calcula comprimento e inclinação de cada elemento e indica se
%   há componente de viga, ou seja, se é um pórtico
data.L = zeros(Nel, 1);
data.Q = zeros(Nel, 1);
for e = 1:Nel
    
    x1 = coor(con(e, 1), 1);
    x2 = coor(con(e, 2), 1);  
    y1 = coor(con(e, 1), 2);
    y2 = coor(con(e, 2), 2);
    
    data.L(e) = sqrt((x2-x1)^2 + (y2-y1)^2);
    data.Q(e) = atan2(y2-y1, x2-x1)*180/pi;
    
    if ismember(e, trelicas)
        data.viga(e) = 0;
    else
        data.viga(e) = 1;
    end
        
end


% Função para cálculo da matriz de rigidez
Ke = @(E, A, I, L, viga)...
    [ ...
     1  0  0 -1  0  0  ;
     0  0  0  0  0  0  ;
     0  0  0  0  0  0  ;
    -1  0  0  1  0  0  ;
     0  0  0  0  0  0  ;
     0  0  0  0  0  0  ;
    ] * E*A/L + ...
    [ ...
     0      0      0      0     0     0      ;
     0      12     6*L    0    -12    6*L    ;
     0      6*L    4*L^2  0    -6*L   2*L^2  ;
     0      0      0      0     0     0      ;
     0     -12    -6*L    0     12   -6*L    ;
     0      6*L    2*L^2  0    -6*L   4*L^2  ;
    ] * E*I/(L^3) * viga;

% Função para cálculo da matriz de massa
Me = @(rho, A, L, viga)...
    [ ...
     2  0  0  1  0  0  ;
     0  0  0  0  0  0  ;
     0  0  0  0  0  0  ; 
     1  0  0  2  0  0  ;
     0  0  0  0  0  0  ;
     0  0  0  0  0  0  ;
    ] * rho*A*L/6 + ...
    [ ...
     0       0       0       0      0      0      ;
     0       156     22*L    0      54    -13*L   ;
     0       22*L    4*L^2   0      13*L  -3*L^2  ;
     0       0       0       0      0      0      ;
     0       54      13*L    0      156   -22*L   ;
     0      -13*L   -3*L^2   0     -22*L   4*L^2  ;
    ] * rho*A*L/420 * viga;

% Função para cálculo da matriz de transformação
Te = @(theta)...
    [
     cosd(theta)  sind(theta)  0          0           0            0  ;
    -sind(theta)  cosd(theta)  0          0           0            0  ;
     0            0            1          0           0            0  ;
     0            0            0          cosd(theta) sind(theta)  0  ;
     0            0            0         -sind(theta) cosd(theta)  0  ;
     0            0            0          0            0           1  ;
    ];

% Função para cálculo do vetor de carregamento distribuído
Fe = @(q, L) q/12*[0 6*L L^2 0 6*L -L^2]';

% Itera cada elemento e monta matrizes globais de massa e rigidez e vetor
%   de carregamento
Kg = zeros(Ngdl);
Mg = zeros(Ngdl);
Fg = zeros(Ngdl, 1);
for e = 1:Nel
    
    ke = Ke(data.E(e), data.A(e), data.I(e), data.L(e), data.viga(e));
    me = Me(data.rho(e), data.A(e), data.L(e), data.viga(e));
    
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
    
    if (nod2 <= pontos.C.nod) ||( ...
        (nod1 >= pontos.E.nod) && (nod2 <= pontos.H.nod))
        
        fe  = Fe(-80*9.8/9, data.L(e));
        Fg(3*nod1-2:3*nod1) = Fg(3*nod1-2:3*nod1) + fe(1:3);
        Fg(3*nod2-2:3*nod2) = Fg(3*nod2-2:3*nod2) + fe(4:6);
        
    end
    
end

% Cria matrizes de rigidez e massa apenas com os G.L. livres
list = 1:Ngdl;
id_fix = [3*pontos.I.nod-2 3*pontos.I.nod-1 3*pontos.G.nod-2:3*pontos.G.nod];
id_free = list(ismember(list, id_fix) == 0);
Kgm = Kg(id_free, id_free);
Mgm = Mg(id_free, id_free);

end