function plot_struct(coor, con)

Nel = size(con, 1);

for e = 1:Nel
    
    x1 = coor(con(e, 1), 1);
    x2 = coor(con(e, 2), 1);
    
    y1 = coor(con(e, 1), 2);
    y2 = coor(con(e, 2), 2);
    
    if y1 == y2
        cor = 'b';
    elseif x1 == x2
        cor = 'r';
    else
        cor = 'g';
    end
    plot([x1 x2], [y1 y2], strcat('.-',cor)); hold on;
    
end
axis equal;