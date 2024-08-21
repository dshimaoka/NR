function d = angdiff(th1, th2)
  
    switch nargin
        case 1
            if length(th1) == 2
                d = th1(1) - th1(2);
            else
                d = th1;
            end
        case 2
            if length(th1) > 1 && length(th2) > 1
                % if both arguments are vectors, they must be the same
                assert(all(size(th1) == size(th2)), 'SMTB:angdiff:badarg', 'vectors must be same shape');
            end
            % th1 or th2 could be scalar
            d = th1 - th2;
    end
    
    % wrap the result into the interval [-pi pi)
    d = mod(d+pi, 2*pi) - pi;
end