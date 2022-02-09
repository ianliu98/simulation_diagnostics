function fre = Hikishima_frequency(by, bz, dt)

len = length(by);
fre = ones(len,1) * 1e-30;

pi2 = 2.0 * pi;
pi5_1 = pi/6.0;
pi5_0 = pi2 - pi5_1;
dt_inv = 1.0 / dt;

rad = atan2(bz, by);
rad = (rad <= 0) * pi2 + rad;
for i=2:len
    radi = rad(i);
    
    if (rad(i-1) >= pi5_0)
        if ((radi <= pi5_1) && (radi >= 0))
            radi = radi + pi2;
        end
    end
    
    if ((radi >= pi5_0) && (rad(i-1) <= pi5_1))
        radi = radi - pi2;
    end
    
    fre(i) = (radi - rad(i-1)) * dt_inv;
end

fre(1) = NaN;

end