% test on |u_1| , where u ~ Unif(S^{n-1})

NVEC = ceil(2.^(5:0.2:11));
u1vec = zeros(size(NVEC));
j = 1;
for n = NVEC
    nexp = 100*n;
    u1abs = 0;
    for i = 1:nexp
        u = randn(n,1);
        u = u/norm(u);
        u1abs = u1abs + abs(u(1));
    end
    u1abs = u1abs/nexp;
    u1vec(j) = u1abs;
    disp(j);
    j = j+1;
    
end

plot(NVEC, u1vec.*sqrt(NVEC));