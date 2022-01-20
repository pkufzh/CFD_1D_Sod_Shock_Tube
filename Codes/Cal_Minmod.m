
% Function Module: Cal. minmod(a, b)
% Or minmod(a, b) = 0.5 * (sign(a) + siagn(b)) * min(abs(a), abs(b))
% a, b is 1-D array with same length

function Q = Cal_Minmod(a, b)
    
    Np = length(a);
    Q = zeros(1, Np);

    for i = 1 : Np
        if ((a(i) * b(i)) > 0)
            if (abs(a(i)) > abs(b(i)))
                Q(i) = b(i);
            else
                Q(i) = a(i);
            end
        else
            Q(i) = 0;
        end
    end


end