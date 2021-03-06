function plot_struct(coor, con, str)

% Plota cada elemento da estrutura da ponte
Nel = size(con, 1);
for e = 1:Nel
    
    x1 = coor(con(e, 1), 1);
    x2 = coor(con(e, 2), 1);
    y1 = coor(con(e, 1), 2);
    y2 = coor(con(e, 2), 2);
    plot([x1 x2], [y1 y2], str);
    
end
axis equal;
axis([-5 70 -10 25]);