function hl = plot_ellipse(ax, x0, y0, a,b,  orientation)

majaxis='x';
% in the ellipse coordinate:
if strcmp( majaxis, 'x')
    x = linspace(-a,a, 50);
    y_prime = sqrt(b^2 - b^2/a^2 .* (x ).^2);
    
else strcmp(majaxis,'y')
    x = linspace(-b,b, 50);
    y_prime = sqrt(a^2 - (a/b)^2.*(x).^2);
end
    

y_pos =  y_prime; y_neg = -y_prime;
if ~isrow(y_pos)
    y_pos = y_pos'; y_neg = y_neg';
end
ys = [y_pos, y_neg];

R= rotz(-orientation);
coord_cart_pos = inv(R)*[x; y_pos] + [x0;y0] ;
coord_cart_neg = inv(R)*[x; y_neg] + [x0;y0];

hl(1)= plot(ax, coord_cart_pos(1,:), coord_cart_pos(2,:),'-k','linewidth',1.2);
hl(2)= plot(ax, coord_cart_neg(1,:), coord_cart_neg(2,:),'-k','linewidth',1.2);

%hl(2) = plot(ax, x, y_neg, '-k', 'linewidth',1.2);




return