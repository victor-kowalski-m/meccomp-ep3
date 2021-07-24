clear
clc
close all

dx = 1;

pA.coor = [0  0 ]; pA.con = 1; % A
pB.coor = [5  0 ]; % B
pC.coor = [9  0 ]; % C
pD.coor = [18 0 ]; % D
pE.coor = [36 0 ]; % E
pF.coor = [36 18]; % F
pG.coor = [36 -5]; % G
pH.coor = [65  0]; % H
 
idx = 1;
coor = [pA.coor];
con = [];

coor_viga_hor = calcula_pontos_internos(pA.coor, pH.coor, dx);
Ncoor_viga_hor = size(coor_viga_hor, 1);
for i=2:Ncoor_viga_hor
    coor = [coor coor_viga_hor(i,:)];
    con = [con; [i-1 i]];
    for ponto=[pB pC pD pE pH]
        if isequal(coor_viga_hor(i,:), ponto.coor)
            ponto.con = i;
        end
    end
end

pG.con = i+1;
coor = [coor pG.coor];
coor_viga_ver = calcula_pontos_internos(pG.coor, pF.coor, dx);
Ncoor_viga_ver = size(coor_viga_ver, 1);
pulouE = 0;
for i=2:Ncoor_viga_ver
    if isequal(coor_viga_ver(i,:), pE.coor)
        con = [con; [con(end, 2) pE.con]];
        pulouE = 1;
    else
        coor = [coor coor_viga_ver(i,:)];
        con = [con; [con(end, 2) Ncoor_viga_hor+i-pulouE]];
    end
end

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

