function [autovec, autoval] = analise_modal(coor, con, Kgm, Mgm, id_free, Ngdl, dx)

scale = 10; % define fator de escala das deformaçôes
figure; % cria figura para os plots
 
% Encontra frequências de ressonância e modos de vibrar
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
    title(mod+"º modo")
    hold off

end

end
