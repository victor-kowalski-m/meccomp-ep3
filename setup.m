clear
clc
close all

dx = 3;

pontos.A.coor = [0  0 ]; pontos.A.con = 1; % A
pontos.B.coor = [5  0 ]; pontos.B.con = 0; % B
pontos.C.coor = [9  0 ]; pontos.C.con = 0; % C
pontos.D.coor = [18 0 ]; pontos.D.con = 0; % D
pontos.E.coor = [36 0 ]; pontos.E.con = 0; % E
pontos.F.coor = [36 18]; pontos.F.con = 0; % F
pontos.G.coor = [36 -5]; pontos.G.con = 0; % G
pontos.H.coor = [65  0]; pontos.H.con = 0; % H

nos = 0;

[AB, nAB] = calcula_pontos_internos(pontos.A.coor, pontos.B.coor, dx);
coor = AB;
con = [(1:nAB-1)', (2:nAB)'];
nos = nos + nAB;
pontos.B.con = nos;

ps = ['B', 'C', 'D', 'E', 'H'];
for i=1:length(ps)-1
    [pIn, NpIn] = calcula_pontos_internos( ...
        pontos.(ps(i)).coor, pontos.(ps(i+1)).coor, dx);
    pIn = pIn(2:end, :); NpIn = NpIn - 1;
    coor = [coor; pIn];
    con = [con; nos+[(0:NpIn-1)', (1:NpIn)']];
    nos = nos + NpIn;
    pontos.(ps(i+1)).con = nos;
end

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

for p = ['B', 'D', 'H']
    con(end+1, :) = [pontos.(p).con pontos.F.con];
end


% [pontos, Npontos] = calcula_pontos_internos( ...
%     pG.coor, pontos_hor(i+1).coor, dx);
% pontos = pontos(2:end, :); Npontos = Npontos - 1;
% coor = [coor; pontos];
% con = [con; (nos-1)+[(1:Npontos-1)', (2:Npontos)']];
% nos = nos + Npontos;

% [pontos, Npontos] = calcula_pontos_internos(pB.coor, pC.coor, dx);
% pontos = pontos(2:end, :); Npontos = Npontos - 1;
% coor = [coor; pontos];
% con = [con; (nos-1)+[(1:Npontos-1)', (2:Npontos)']];
% nos = nos + Npontos;
% 
% [CD, Npontos] = calcula_pontos_internos(pB.coor, pC.coor, dx);
% pontos = pontos(2:end, :); Npontos = Npontos - 1;
% coor = [coor; CD];
% con = [con; (nos-1)+[(1:Npontos-1)', (2:Npontos)']];
% nos = nos + Npontos;


% coor_viga_hor = calcula_pontos_internos(pA.coor, pH.coor, dx);
% Ncoor_viga_hor = size(coor_viga_hor, 1);
% for i=2:Ncoor_viga_hor
%     coor = [coor coor_viga_hor(i,:)];
%     con = [con; [i-1 i]];
%     for ponto=[pB pC pD pE pH]
%         if isequal(coor_viga_hor(i,:), ponto.coor)
%             ponto.con = i;
%         end
%     end
% end
% 
% pG.con = i+1;
% coor = [coor pG.coor];
% coor_viga_ver = calcula_pontos_internos(pG.coor, pF.coor, dx);
% Ncoor_viga_ver = size(coor_viga_ver, 1);
% pulouE = 0;
% for i=2:Ncoor_viga_ver
%     if isequal(coor_viga_ver(i,:), pE.coor)
%         con = [con; [con(end, 2) pE.con]];
%         pulouE = 1;
%     else
%         coor = [coor coor_viga_ver(i,:)];
%         con = [con; [con(end, 2) Ncoor_viga_hor+i-pulouE]];
%     end
% end

% 
% con_viga_hor = (1:Ncoor_viga_hor) + (idx-1);
% coor = [coor; coor_viga_hor];
% con =  [con; [con_viga_hor(1:end-1)', con_viga_hor(2:end)']];
% elems_viga_hor = [idx idx+Ncoor_viga_hor];
% idx = elems_viga_hor(2);
% 
% coor_viga_ver = calcula_pontos_internos(pontos.G, pontos.F, dx);
% Ncoor_viga_ver = size(coor_viga_ver, 1);
% con_viga_ver = (1:Ncoor_viga_ver) + (idx-1);
% coor = [coor; coor_viga_ver];
% con =  [con; [con_viga_ver(1:end-1)', con_viga_ver(2:end)']];
% elems_viga_ver = [idx idx+Ncoor_viga_ver];
% idx = elems_viga_ver(2);


% pontos_desenho = [ 
%     0  0 ; % A
%     5  0 ; % B
%     9  0 ; % C
%     18 0 ; % D
%     36 0 ; % E
%     36 18; % F
%     36 -5; % G
%     65  0; % H
% ];

% con = [
%     
%     % Viga horizontal
%     1  8  1;
%     
%     % Viga vertical
%     7  6 1;
%     
%     % Treliças
%     2  6  0;
%     4  6  0;
%     6  8  0;
% ];

%[coor, con] = discretiza(coorDesenho, conDesenho);
plot_struct(coor, con)

