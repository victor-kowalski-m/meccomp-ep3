function transiente...
    (coor, con, pontos, data, Kgm, Mgm, Fg, id_free, Ngdl)

figure('Name','Análise Transiente','NumberTitle','off');

Fgm = Fg(id_free); % Vetor de força nos GL nao nulos

% Rayleigh
a0 = 0.1217;
a1 = 0.0087;
Cgm = a0*Mgm + a1*Kgm;

% Método da aceleração média constante
gamma = 0.5;
beta = 0.25;

% Roda a análise transiente para cada valor de va e dt
vas = [0.5 1 2];
dts = {0.01 [0.1 0.01 0.005] 0.01};
for j=1:length(vas)
    for dt=dts{j}

    % Define va e parâmetros de tempo da iteração
    va = vas(j);
    TF = 55/va;
    t = 0:dt:TF;

    % Define matrizes/vetores do método Newmark Beta
    Meq = Mgm + dt*gamma*Cgm + dt^2*beta*Kgm;
    U1 = Kgm\(20*Fgm);
    V1 = zeros(length(U1), 1);
    A1 = zeros(length(U1), 1);
    
    V = zeros(Ngdl, 1); % vetor de cargas verticais Vi
    
    % Vetores auxiliares para a plotagem
    U_plot = zeros(Ngdl, 1);
    deslocA = zeros(length(t), 1);
    deslocB = zeros(length(t), 1);
    deslocC = zeros(length(t), 1);
    deslocF = zeros(length(t), 1);

    for i = 1:length(t)
        
        F2 = Fg; % reseta vetor de forças para vetor genérico

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
        F2(3*pontos.E.nod-2:3*pontos.H.nod) = ...
            F2(3*pontos.E.nod-2:3*pontos.H.nod)*N2;

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

        % Cálcula grandezas pelo método Newmark Beta
        Feq = F2(id_free) - Cgm*(V1 + dt*(1-gamma)*A1) ...
            - Kgm*(U1 + dt*V1 + dt^2/2*(1-2*beta)*A1);
        A2 = Meq\Feq;
        U2 = U1 + dt*V1 + dt^2/2*((1-2*beta)*A1 + 2*beta*A2);
        V2 = V1 + dt*((1-gamma)*A1 + gamma*A2);

        % Armazena deslocamentos da estrutura nos pontos determinados
        U_plot(id_free) = U2;
        deslocA(i) = U_plot(3*pontos.A.nod-1);
        deslocB(i) = U_plot(3*pontos.B.nod-1);
        deslocC(i) = U_plot(3*pontos.C.nod-1);
        deslocF(i) = U_plot(3*pontos.F.nod-2);

        % Atualiza valores da iteração anterior 
        U1 = U2;
        V1 = V2;
        A1 = A2;

    end

    % Se va for 1, plota deslocamento do ponto A para diferentes dts
    if va == 1
        subplot(2, 2, 1);
        plot(t, deslocA)
        hold on
        legend(string(dts{j})+"s")
        title('\bf{Ponto A} ({\boldmath$v_a='+string(va)+'$})', ...
            'Interpreter','latex')
        ylabel("U (m)")
        xlabel("t (s)")
    end

end

% Plota deslocamento dos pontos determinados para cada va e dt = 0.01
subplot(2, 2, j+1);
hold on
plot(t, deslocA)
plot(t, deslocB)
plot(t, deslocC)
plot(t, deslocF)
hold off
legend(["A" "B" "C" "F"])
title('\bf{Pontos ABCF} ({\boldmath$v_a='+string(va)+'$})', ...
    'Interpreter','latex')
ylabel("U (m)")
xlabel("t (s)")
end

end