function dfdx = four_order_appr(f,dx)
    n = length(f);
    dfdx = zeros(n-4,1);
    for i=3:n-2
        dfdx(i-2) = ( f(i-2) - 8.0 * f(i-1) + 8.0 * f(i+1) - f(i+2) ) / (12.0 * dx);
    end
end