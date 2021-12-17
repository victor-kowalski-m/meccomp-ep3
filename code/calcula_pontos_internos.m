function [pontos_in, Npontos] = calcula_pontos_internos(ponto1, ponto2, dx)
    
    % Cálcula distância entre os dois pontos e número de pontos entre eles
    vetor = ponto2 - ponto1;
    norma = norm(vetor);
    Npontos = fix(norma/dx)+1;
    pontos_in = zeros(Npontos, 2);
    
    % Adiciona pontos entre os dois iniciais, conforme o dx
    for i=1:Npontos
        pontos_in(i, :) = ponto1 + (i-1)*dx*vetor/norma;
    end
    
    % Se sobrou espaço menor que dx ao final, utiliza um passo menor
    dist_ponto_final = norm(ponto2 - pontos_in(end, :));
    if dist_ponto_final
        pontos_in(end+1, :) = pontos_in(end, :) + ...
            dist_ponto_final*vetor/norma;
        Npontos = Npontos + 1;
    end
    
end