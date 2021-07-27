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
pontos.H.coor = [45  0]; pontos.H.nod = 0; % H
pontos.I.coor = [65  0]; pontos.I.nod = 0; % I

coor = pontos.A.coor; con = []; nos = 1; % inicialização de matrizes de
    % coordenadas e conectividade

%%% Insere viga horizontal nas matrizes de coord. e conec.
ps = ['A', 'B', 'C', 'D', 'E', 'H', 'I'];
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
for p = ['B', 'D', 'I']
    con(end+1, :) = [pontos.(p).nod pontos.F.nod];
end
trelicas = (vigaVer(end)+1):size(con, 1);

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

data.I = zeros(Nel, 1);
data.I(vigaHor) = 0.9*0.9^3/12;
data.I(vigaVer) = 1.8^3*0.9/12;
data.I(trelicas) = pi*0.05^4/32;

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


%%% Funções para cálculo das matrizes de rigidez, massa e transformação
Ke = @(E, A, I, L, viga)...
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
    ] * E*I/(L^3) * viga;

Me = @(rho, A, L, viga)...
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
    ] * rho*A*L/420 * viga;

Te = @(theta)...
    [
     cosd(theta)  sind(theta)  0          0           0            0;
    -sind(theta)  cosd(theta)  0          0           0            0;
     0            0            1          0           0            0;
     0            0            0          cosd(theta) sind(theta)  0;
     0            0            0         -sind(theta) cosd(theta)  0;
     0            0            0          0            0           1
    ];

Fe = @(q, L) q/12*[0 6*L L^2 0 6*L -L^2]';

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
   
    if nod2 <= pontos.C.nod ||( ...
            (nod1 >= pontos.E.nod) && (nod2 <= pontos.H.nod))
        
        fe  = Fe(-80*9.8/9, data.L(e));
        Fg(3*nod1-2:3*nod1) = Fg(3*nod1-2:3*nod1) + fe(1:3);
        Fg(3*nod2-2:3*nod2) = Fg(3*nod2-2:3*nod2) + fe(4:6);
        
    end
    
end

list = 1:Ngdl;
id_fix = [3*pontos.I.nod-2 3*pontos.I.nod-1 3*pontos.G.nod-2:3*pontos.G.nod];
id_free = list(ismember(list, id_fix) == 0);
Ngdla = length(id_free);

Kgm = Kg(id_free, id_free);
Mgm = Mg(id_free, id_free);
Fgm = Fg(id_free);

a0 = 0.1217;
a1 = 0.0087;
Cgm = a0*Mgm + a1*Kgm;
gamma = 0.5;
beta = 0.25;

U_plot = zeros(Ngdl, 1);
scale = 500;

vas = [0.5 1 2];
dts{1} = 0.01;
dts{2} = [0.1 0.01 0.005];
dts{3} = 0.01;

for j=1:length(vas)%va=[0.5 1 2]
for dt=dts{j}
    
va = vas(j);
TF = 55/va;
t = 0:dt:TF;

Meq = Mgm + dt*gamma*Cgm + dt^2*beta*Kgm;

U1 = Kgm\(20*Fgm);
V1 = zeros(Ngdla, 1);
A1 = zeros(Ngdla, 1);
V = zeros(Ngdl, 1);

deslocA = zeros(length(t), 1);
deslocB = zeros(length(t), 1);
deslocC = zeros(length(t), 1);
deslocF = zeros(length(t), 1);

for i = 1:length(t)
        
    F2 = Fg;

    % Define N1
    if t(i) > 27/va
        if t(i) <= 47/va
            N1 = va*t(i)-7;
        else
            N1 = 40;
        end
    else
        N1 = 20;
    end
    F2(1:3*pontos.C.nod) = F2(1:3*pontos.C.nod)*N1;
    
    % Define N2
    if t(i) <= 20/va
        N2 = 20-va*t(i);
    else
        N2 = 0;
    end
    F2(3*pontos.E.nod-2:3*pontos.H.nod) = F2(3*pontos.E.nod-2:3*pontos.H.nod)*N2;
    
    % Define Vi
    for e=(pontos.C.nod+1):(pontos.E.nod-1)
        if t(i) > (pontos.E.coor(1)-coor(e, 1))/va
            if t(i) < (pontos.E.coor(1)+20-coor(e, 1))/va
                V(3*e-1) = -80*9.8*(1-cos(2*pi*va*t(i)))/2;
            else
                V(3*e-1) = 0;
            end
        else
            V(3*e-1) = 0;
        end
    end
    F2 = F2 + V;
    
%     F3 = F2;
    F2 = F2(id_free);
    
    Feq = F2 - Cgm*(V1 + dt*(1-gamma)*A1) - Kgm*(U1 + dt*V1 + dt^2/2*(1-2*beta)*A1);
    A2 = Meq\Feq;
    U2 = U1 + dt*V1 + dt^2/2*((1-2*beta)*A1 + 2*beta*A2);
    V2 = V1 + dt*((1-gamma)*A1 + gamma*A2);
    
    U_plot(id_free) = U2;
    deslocA(i) = U_plot(3*pontos.A.nod-1);
    deslocB(i) = U_plot(3*pontos.B.nod-1);
    deslocC(i) = U_plot(3*pontos.C.nod-1);
    deslocF(i) = U_plot(3*pontos.F.nod-2);
    
%     if mod(i, 20)
%     
%     coorExag = coor + scale*[U_plot(1:3:end) U_plot(2:3:end)];
% %     figure(1);
%     plot_struct(coorExag, con, '-r');
%     axis([-10 80 -10 25])
%     pause(0.001);
%     grid;
%     clf;
% %     figure(2);
% %     plot(F3(2:3:end));
% %     axis([-50 300 -3000 3000])
% %     pause(0.001);
%     end
    
    U1 = U2;
    V1 = V2;
    A1 = A2;
       
end

if va == 1
    subplot(2, 2, 1);
    plot(t, deslocA)
    hold on
    legend(string(dts{j})+"s")
    title('\bf{Ponto A} ({\boldmath$v_a='+string(va)+'$})','Interpreter','latex')
end

end

subplot(2, 2, j+1);
hold on
plot(t, deslocA)
plot(t, deslocB)
plot(t, deslocC)
plot(t, deslocF)
hold off
legend(["A" "B" "C" "F"])
title('\bf{Pontos ABCF} ({\boldmath$v_a='+string(va)+'$})','Interpreter','latex')

end
