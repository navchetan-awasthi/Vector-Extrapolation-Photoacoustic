% Extrapolate iterative method F using specified method.
% Use L cycles of k-order extrapolation, starting from x0.
function [x0, residuals,final1] = my_extrapolate_tv_hair(V,S,A_b,b,x0, k, L, method,tol)
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
    for t = 1:L
        Q(:, 1)=tv(b,x0,A_b,V,S);   
        data_diff=(b-A_b*x0);
        pj_error=[pj_error sum(abs(data_diff.^2))];
        data_diff=(b-A_b*Q(:,1));
        pj_error=[pj_error sum(abs(data_diff.^2))];
        residuals(t)=(pj_error(end-1)-pj_error(end))*100/pj_error(end-1);
        residuals(t)';
        if abs(residuals(t))<tol
            break
        end
        for i = 1:k 
            Q(:, i+1)=tv(b,Q(:,i),A_b,V,S);   
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


% Total Variation Method
function [x] = tv(y,x,A_b,V,S)
    mu = 1e-3;
    matA=A_b;
    tau=1e-2;
    [M,N] = size(A_b);
    inverse_term=0;
    invLS = @(x) inverse_term*x;        
    phi = @(x) TVnorm(x);
    P = @(x) x;
    PT = @(x) x;
    PTx = PT(x);
    u = PTx;
    bu = 0*u;
    threshold = tau/mu;
    criterion(1) = 1;
    % Compute and store initial value of the objective function
    numA = 0;
    A=A_b;
    resid =  y-A*x;
    numA = numA + 1;
    prev_f = 0.5*(resid(:)'*resid(:)) + tau*phi(reshape(u,200,200));
    times(1) = 0;
    objective(1) = prev_f;
    puy=zeros(200,200);
    pux=puy;
    maxiter=1;
    isTVinitialization=1;
    TViters=10;
    stopCriterion=2;
    tolA=1e-3;
    for outer = 1:maxiter
    PTx=x;
    [u,pux,puy] = chambolle_prox_TV_stop(real(PTx-bu), 'lambda', threshold, 'maxiter', TViters, 'dualvars',[pux puy]);
    r = A_b'*y + mu*(reshape(u,40000,1)+bu);
    size(r);
    sigma = diag(S);
    sigma(10:40000)=0;
    V1 = V;
    alpha=0.01;
    nume = 1;
    deno = sigma.^2 + alpha^2;
    b = r;
    beta = V'*b;
    phi = nume./deno;
    zeta = phi.*beta;
    x = V*(zeta);
    end   
     
end



