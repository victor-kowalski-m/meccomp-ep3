function [autovec, autoval] = modal...
    (coor, con, Kgm, Mgm, id_free, Ngdl, dx)

scale = 10; % define fator de escala das deformaçôes
figure('Name',"Análise Modal (dx="+dx+")",'NumberTitle','off'); 
 
% Encontra as 6 primeiras frequências de ressonância e modos de vibrar
[autovec, autoval] = eigs(Mgm\Kgm, 6, 'SM');
autoval = sqrt(diag(autoval))/(2*pi);

U = zeros(Ngdl, 1); % vetor de deslocamentos 

% Plota cada modo de vibrar da estrutura
for mod=1:6

    freq = autoval(mod);
    U(id_free) = autovec(:, mod);
    coorExag = coor + scale*[U(1:3:end) U(2:3:end)];

    subplot(2, 3, mod)
    hold on
    plot_struct(coorExag, con, '-r');
    plot_struct(coor, con, '-b');
    title(mod+"º modo ("+autoval(mod)+" Hz)")
    hold off

end

end
