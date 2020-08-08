% Extrapolate iterative method F using specified method.
% Use L cycles of k-order extrapolation, starting from x0.
function [x0, residuals,final1] = my_extrapolate_rsd(A_b,b,x0, k, L, method,tol)
tic
    N = numel(x0);
    Q = zeros(N, k+1);
    switch upper(method)
        case 'MPE', method = @mpe;
        case 'RRE', method = @rre;
        case 'N/A', method = '';
        case '', warning('Extrapolate:PlainIteration', 'No extrapolation');
        otherwise, error('Extrapolate:UnknownMethod', method);
    end
    % Perform L cycles of extrapolation method
        alpha=1;
        q=4;
    residuals = zeros(L, 1);
    pj_error=[];
            data_diff=(b-A_b*x0);
        pj_error=[pj_error sum(abs(data_diff.^2))];
        r0=A_b*x0-b;
        l=A_b'*r0+alpha*x0;
        fac=(norm(l,2)^2)/(norm(A_b*l,2)^2+alpha*norm(l,2)^2);
        Q(:, 1)=x0;       
    for t = 1:L       
        for i = 1:k 
            alpha=alpha/q;
            r0=A_b*Q(:, i)-b;
            l=A_b'*r0+alpha*x0;
            fac=(norm(l,2)^2)/(norm(A_b*l,2)^2+alpha*norm(l,2)^2);
            Q(:, i+1)=Q(:, i)-fac*l;
        end       
        if isempty(method) % No extrapolation.
            x0 = Q(:, end); % Just take the last vector.
            continue;
        end
        % Compute differences (k+1)
        for i = k:-1:1 
            Q(:, i+1) = Q(:, i+1) - Q(:, i);
        end
        Q(:, 1) = Q(:, 1) - x0;
        % Perform QR decomposition
        [Q, R] = MGS(Q); 
        % Perform extrapolation
        [gamma] = method(R, k); % s.t. x0 = X * gamma
        xi = 1 - cumsum(gamma(1:k)); % s.t. x0' = x0 + U * xi
        eta = R(1:k, 1:k) * xi(:); % since U = Q * R
        x0 = x0 + Q(:, 1:k) * eta; % s.t. x0' = x0 + Q * R * xi   
        data_diff=(b-A_b*x0);
        pj_error=[pj_error sum(abs(data_diff.^2))];
        residuals(t)=(pj_error(end-1)-pj_error(end))*100/pj_error(end-1);
        residuals(t)'
        if abs(residuals(t))<tol
            break
        end
        residuals(t);
        Q(:, 1)=x0;
    end
    final1=toc;
end

% Minimal Polynomial Extrapolation
function [gamma, residual] = mpe(R, k)
    c = backsub(R(1:k, 1:k), -R(1:k, k+1));
    c = [c; 1];
    gamma = c / sum(c);
    residual = abs(gamma(end)) * R(end, end);
end

% Reduced Rank Extrapolation
function [gamma, residual] = rre(R, k)
    e = ones(k+1, 1);
    d = backsub(R, backsub(R', e));  
    lambda = 1 / sum(d);
    gamma = lambda * d;
    residual = sqrt(lambda);
end

