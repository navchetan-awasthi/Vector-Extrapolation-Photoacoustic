% Extrapolate iterative method F using specified method.
% Use L cycles of k-order extrapolation, starting from x0.
function [x0, residuals,final1] = my_extrapolate_tv(A_b,b,inverse_term,x0, k, L, method,tol)
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
    residuals = zeros(L, 1);
    pj_error=[];
            data_diff=(b-A_b*x0);
        pj_error=[pj_error sum(abs(data_diff.^2))];
    for t = 1:L
        t;
        % Compute (k+1) vectors, in addition to x0
        Q(:, 1)=tv(b,x0,A_b,inverse_term);   
%         residuals(t) = norm(A_b*x0 - A_b*Q(:, 1), 2)/norm(A_b*Q(:,1)-b,2);

        data_diff=(b-A_b*Q(:,1));
        pj_error=[pj_error sum(abs(data_diff.^2))];
        residuals(t)=(pj_error(end-1)-pj_error(end))*100/pj_error(end-1);
        residuals(t)';
        if abs(residuals(t))<tol
            break
        end
        for i = 1:k 
            Q(:, i+1)=tv(b,Q(:,i),A_b,inverse_term);  
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
        residuals(t)';
        if abs(residuals(t))<tol
            break
        end
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


% total variation
function [x] = tv(y,x,A_b,inverse_term)
    mu = 1e-3;
    tau=1e-2;
    invLS = @(x) inverse_term*x;        
    PT = @(x) x;
    PTx = PT(x);
    u = PTx;
    bu = 0*u;
    threshold = tau/mu;
    puy=zeros(201,201);
    pux=puy;
    maxiter=2;
    TViters=10;
    for outer = 1:maxiter
    PTx=x;
    [u,pux,puy] = chambolle_prox_TV_stop(real(PTx-bu), 'lambda', threshold, 'maxiter', TViters, 'dualvars',[pux puy]);
    r = A_b'*y + mu*(reshape(u,40401,1)+bu);
    x = invLS(r);
    end   
     
end



