function [pontos_in, Npontos] = calcula_pontos_internos(ponto1, ponto2, dx)
    
    vetor = ponto2 - ponto1;
    norma = norm(vetor);
    Npontos = fix(norma/dx)+1;
    pontos_in = zeros(Npontos, 2);
    
    for i=1:Npontos
        pontos_in(i, :) = ponto1 + (i-1)*dx*vetor/norma;
    end
    
    dist_ponto_final = norm(ponto2 - pontos_in(end, :));
    if dist_ponto_final
        pontos_in(end+1, :) = pontos_in(end, :) + dist_ponto_final*vetor/norma;
        Npontos = Npontos + 1;
    end
    
end