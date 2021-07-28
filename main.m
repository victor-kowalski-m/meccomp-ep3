clear
clc
close all

for dx=[3 2 1]
    
    % Monta a estrutura e matrizes de MEF
    [coor, con, pontos, data, Kgm, Mgm, Fg, id_free, Ngdl] = setup(dx);
    
    % Análise modal para cada dx
    [autovec, autoval] = modal(coor, con, Kgm, Mgm, id_free, Ngdl, dx);
    fprintf("As 6 primeiras frequências naturais para dx = %d \n", dx);
    disp(autoval);
    
end % estrutura usada nos próximos será com dx=1

% Análise transiente
transiente(coor, con, pontos, data, Kgm, Mgm, Fg, id_free, Ngdl);

% % Análise harmônica
harmonica(coor, con, pontos, data, Kgm, Mgm, Fg, id_free, Ngdl);