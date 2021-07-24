function [coorDisc, conDisc] = discretiza(coor, con, dx)
    Nel = size(con, 1);
    for i=1:Nel
        if con(i, 3)
            pontos = calcula_pontos(coor, con(i, 1), con(i, 2), dx);
        end
    end
end