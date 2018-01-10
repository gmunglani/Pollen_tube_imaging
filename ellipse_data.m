function [tip_final,center,phi,axes,tip_check,n,flag_tol] = ellipse_data(tip_new)

    tip_news = zeros(length(tip_new),2);
    tip_news(:,2) = smooth(tip_new(:,2));
    tip_news(:,1) = smooth(tip_new(:,1));
        
    x = tip_news(:,1);
    y = tip_news(:,2);
    
    a = ellipse_fit(x,y);
    center = ellipse_center(a);
    axes = ellipse_axis_length(a); 
    [phi n] = ellipse_angle_of_rotation2(a,axes);
   
    flag_tol = isreal(axes);
    if (flag_tol == 1)  
        [xe1, ye1, xe2, ye2] = ellipse_tip(center, axes, phi);
    
        diste1 = pdist2([x y],[xe1 ye1],'euclidean');
        diste2 = pdist2([x y],[xe2 ye2],'euclidean');

        dsum1 = sum(diste1);
        dsum2 = sum(diste2);
    
        if (dsum1 < dsum2)
            tip_check  = [xe1 ye1];
            [pt pos] = min(diste1);
        else
            tip_check  = [xe2 ye2];
            [pt pos] = min(diste2);
        end   
        tip_final = tip_new(pos,:);
    else
        tip_check = [0 0];
        tip_final = [0 0];
    end
end

function a = ellipse_fit(x, y)
    D = horzcat(x .* x, x .* y, y .* y, x, y, ones(size(x,1),1));
    S = D' * D;
    C = zeros(6, 6);
    C(1, 3) = 2;
    C(3, 1) = 2;
    C(2, 2) = -1;
    E = eig(inv(S) * C);
    [V,D] = eig(inv(S) * C);
    [nval n] = max(abs(E));
    a = V(:, n);
end

function cen = ellipse_center(a)
    b = a(2)*0.5;
    c = a(3); 
    d = a(4)*0.5;
    f = a(5)*0.5;
    g = a(6);
    a = a(1);
    num = b * b - a * c;
    x0 = (c * d - b * f) / num;
    y0 = (a * f - b * d) / num;
    cen = [x0 y0];
end
    
function ell = ellipse_axis_length(a)
    b = a(2)*0.5;
    c = a(3); 
    d = a(4)*0.5;
    f = a(5)*0.5;
    g = a(6);
    a = a(1);
    up = 2 * (a * f * f + c * d * d + g * b * b - 2 * b * d * f - a * c * g);
    down1 = (b * b - a * c) * ((c - a) * sqrt(1 + 4 * b * b / ((a - c) * (a - c))) - (c + a));
    down2 = (b * b - a * c) * ((a - c) * sqrt(1 + 4 * b * b / ((a - c) * (a - c))) - (c + a));
    res1 = sqrt(up / down1);
    res2 = sqrt(up / down2);
    ell = [res1 res2];
end

function [ang n] = ellipse_angle_of_rotation2(a, axes)
    b = a(2)*0.5;
    c = a(3); 
    d = a(4)*0.5;
    f = a(5)*0.5;
    g = a(6);
    a = a(1);
    n = 0;
    if (b == 0)
        if (a > c)
            ang = 0;
        else
            ang = pi*0.5;
        end
    else
        if (a > c)
            n = 1;
            ang = atan2(2 * b,(a - c))/ 2;
            if (axes(1) > axes(2)) ang = ang - pi/2; end
        else
            n = 2;
            ang = atan2(2 * b,(a - c)) / 2 - pi/2 ;
            if (axes(1) > axes(2)) ang = ang + pi/2; end
        end
    end
end

function [xe1, ye1, xe2, ye2] = ellipse_tip(center, axes, phi)

if (axes(2) > axes(1))
    a = axes(1);
    b = axes(2);
else
    a = axes(2);
    b = axes(1);
end

cosphi = cos(pi/2);
sinphi = sin(pi/2);

R = [ cos(phi)   -sin(phi)
    sin(phi)   cos(phi)];

xy = [a*cosphi; b*sinphi];
xy = R*xy;

xe1 = xy(1,:) + center(1);
ye1 = xy(2,:) + center(2);

xe2 = -xy(1,:) + center(1);
ye2 = -xy(2,:) + center(2);

end
