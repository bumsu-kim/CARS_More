function Result = NesterovRS(fparam, param)
algname = 'CARS';
%% INITIALIZATION
Result = struct;
n = param.n;
eps = 1e-6; maxit = 100;

x = zeros(n,1); % initial sol (x0)
randAlg = 'U'; % uniform
verbose = false;


if isfield(param, 'eps')
    eps = param.eps;
end
if isfield(param, 'x0')
    x = param.x0;
end
if isfield(param, 'maxit')
    maxit = param.maxit;
end
if isfield(param, 'randAlg')
    randAlg = param.randAlg;
end
if isfield(param, 'verbose')
    verbose = param.verbose;
end
if isfield(fparam, 'fmin')
    fmin = fparam.fmin;
end
f = fparam.f;
objval_seq = zeros(maxit+1,1);
gamma_seq = zeros(maxit+1,1);
objval_seq(1) = f(x); %initialization
num_queries = zeros(maxit+1,1);
num_queries(1) = 1;
damped_cnt = 0;

%% ITERATION
mu_init = 1e-4;
mu_max = 1;
mu = mu_init;
mu_cnt = 0;
gd_cnt = 0;

mu_seq = zeros(maxit+1,1);
mu_seq(1) = mu;
alpha = 1/4/(n+4);0.5;

n_dir_grad = 1; ceil(n/2);
been_there = false;
stuck_count = 0;
CARScounter = zeros(4,1);
for k=1:maxit
    num_queries(k+1) = num_queries(k);
    mu =  1/sqrt(k+1);
    fx = objval_seq(k); % use saved value
    u = PickRandDir(1, n, randAlg)'; 
    u = u/norm(u); % normalize u
    fp = f(x+mu*u);
    num_queries(k+1) = num_queries(k+1) + 1;
    d = (fp - fx)/mu; % directional derivative
    
    delta = -alpha*d*u; % move to the next iterate
    
    fxnew = f(x+delta);
    num_queries(k+1) = num_queries(k+1) + 1; %single query here
    
    fs = [fx, fp, fxnew];
    [fxnew, midx] = min(fs);
    if midx == 1
        % no update
        delta = 0;
    else
        if midx == 2
            delta = mu*u;
        elseif midx == 3
            % delta not changed
        end
    end
    x = x + delta;
%     if norm(delta) < 0.99*mu
%         been_there = true;
%     end
%     if been_there || (fxnew-objval_seq(k) < 1e-6)
%         stuck_count = stuck_count + 1;
%     else
%         mu = mu_init;
%         stuck_count = 0;
%     end
%     if mu>1e-8
%         mu = 0.5* mu + norm(delta) * mu;
%     end
%     mu= 1/sqrt(k+2);
    
    if isnan(fxnew) % should never happen..
        disp('errrrrrr');
    end
%     x = x+delta; % update the solution
    objval_seq(k+1) = fxnew;
    
    if isfield(fparam, 'fmin')
        if (fxnew < fmin + eps) 
            if verbose>1
                disp([algname, ' Converged in ', num2str(k),' steps. Exit the loop']);
                disp(['Function val = ' , num2str(fxnew)]);
            end
            converged = true;
            break;
        end
    end
    if (num_queries(k+1)>param.MAX_QUERIES)
        break;
    end
    mu_seq(k+1) = mu;
end

if (k>=maxit) || (num_queries(k+1)>param.MAX_QUERIES)
    if verbose>1
        disp([algname, ' did not converge in ', num2str(maxit) , ' steps.']);
    end
    converged = false;
end

num_iter = k;
objval_seq = objval_seq(1:num_iter+1);
% sol_seq = sol_seq(1:num_iter+1);
gamma_seq = gamma_seq(1:num_iter+1);
num_queries = num_queries(1:num_iter+1);
mu_seq = mu_seq(1:num_iter+1);

% put into a struct for output
if isfield(param,'save_x')
    if param.save_x
%         Result.sol = sol_seq;
    end
end
Result.objval_seq = objval_seq;
Result.gamma_seq = gamma_seq;
Result.num_iter = num_iter;
Result.num_queries = num_queries;
Result.converged = converged;
Result.damped_cnt = damped_cnt;
Result.mu_cnt = mu_cnt;
Result.gd_cnt = gd_cnt;
Result.mu_seq = mu_seq;
Result.sol = x;
% Result.k_lambda = k_lambda_condition;
if verbose>1
    disp(CARScounter');
end
end
