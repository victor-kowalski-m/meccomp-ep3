clear
clc
close all

for dx=[3 2 1]
    
    % Monta a estrutura e matrizes de MEF
    [coor, con, pontos, data, Kgm, Mgm, Fg, id_free, Ngdl] = setup(dx);
%     [autovec, autoval] = analise_modal(coor, con, Kgm, Mgm, id_free, Ngdl, dx);
    
end

analise_transiente(coor, con, pontos, data, Kgm, Mgm, Fg, id_free, Ngdl);
% analise_harmonica();