
% Gera vetor de carregamento
Fg = Fg*0;
for e=(pontos.C.nod+1):(pontos.E.nod-1)
    Fg(3*e-1) = -80*9.8;
end
Fgm = Fg(id_free);

% Gera vetor de frequências com discretizacao maior (variavel zoom) nas 
%   proximidades das frequências de ressonância
FF = 15;
df = 0.01;
f = [];
ini = 0;
zoom = 100;
for i=1:length(autoval)
    av = autoval(i);
    f = [f ini:df:(av-df) (av-df+df/zoom):df/zoom:(av+df-df/zoom)];
    ini = (av+df);
end
f = [f ini:df:FF];

% Vetores auxiliares para a plotagem
U_plot = zeros(Ngdl, 1);
deslocA = zeros(length(f), 1);
deslocB = zeros(length(f), 1);
deslocC = zeros(length(f), 1);
deslocF = zeros(length(f), 1);

% Resolve a equação dinâmica para cada frequência
for i = 1:length(f)
        
    U2 = (Kgm - (2*pi*f(i))^2*Mgm)\Fgm;
    
    U_plot(id_free) = abs(U2);
    deslocA(i) = U_plot(3*pontos.A.nod-1);
    deslocB(i) = U_plot(3*pontos.B.nod-1);
    deslocC(i) = U_plot(3*pontos.C.nod-1);
    deslocF(i) = U_plot(3*pontos.F.nod-2);
    
       
end

% Plota norma do deslocamento versus frequência nos pontos determinados
hold on
plot(f, deslocA)
plot(f, deslocB)
plot(f, deslocC)
plot(f, deslocF)
hold off
legend(["A" "B" "C" "F"])
title('Resposta em frequência')
ylabel("||U|| (m)")
xlabel("f (Hz)")
grid();
