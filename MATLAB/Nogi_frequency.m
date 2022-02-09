% functions
function fre = Nogi_frequency(by, bz, dt)
    Len = length(by);
    phi_rad = atan2(bz, by);
    phi_continuous = unwrap(phi_rad) / (2 * pi);
    idx = 1;
    idx_prec = 1;
    period = zeros(Len, 1);
    
    for it=2:Len
        tmp_cycle = phi_continuous(it);
        target_cycle = tmp_cycle + 1;
        
        if (phi_continuous(it) - phi_continuous(it-1) >= 0)
            for idx=idx:Len
                if (phi_continuous(idx) > target_cycle)
                    if (idx == 1)
                        break
                    end
                    idx_prec = (target_cycle - phi_continuous(idx-1)) / (phi_continuous(idx) - phi_continuous(idx-1)) ...
                        + idx-1;
                    break
                end
            end
        else
            for idx=idx:-1:1
                if  (phi_continuous(idx) <= target_cycle)
                    idx_prec = (target_cycle - phi_continuous(idx)) / (phi_continuous(idx+1) - phi_continuous(idx)) ...
                        + idx;
                    break
                end
            end
        end
        
        period(it-1) = idx_prec - it;

        if idx == Len
            break 
        end
    end

    period = period * dt;
    fre = 2 * pi ./ period;
end