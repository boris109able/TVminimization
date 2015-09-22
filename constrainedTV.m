function[x] = constrainedTV(A,f)
%% compute min |x|_TV subject to Ax=f
%% split Bregman method
%% refer to "Split Bregman Methods and Frame Based Image Restoration", Jian-Feng Cai, Stanley Osher and Zuowei Shen
%% pp346 top equations
%% delta_b=delta_c=1
%% lambda=mu=1

[m,n]=size(A);
thresh = 0.00001;
u = zeros(n,1);
d = zeros(n-1,1);
b = zeros(n-1,1);
c = zeros(m,1);
lambda = 1;
mu = 1;
B=zeros(n-1,n);
for i=1:(n-1)
    B(i,i)=-1;
    B(i,i+1)=1;
end
while sum(abs(f-A*u))>thresh
    %[u,flag]=pcg((mu*A'*A+lambda*B'*B),(mu*A'*(f-c)+lambda*B'*(d-b)),1e-10,1000000,[],[],u);
    u = (mu*A'*A+lambda*B'*B)\(mu*A'*(f-c)+lambda*B'*(d-b));
    d = wthresh(B*u+b, 's', 1/lambda);
    b = b + B*u-d;
    c = c + A*u-f;
end
x = u;
end