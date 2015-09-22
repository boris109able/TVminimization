%% main function

n=100;

x_star = zeros(n,1);
prob = zeros(n,n);
%load('prob_v3.mat');%continue filling prob from last task

thresh = 0.001;
runNumber = 100;
for K=0:n-1
    fprintf('Computing sparisity level %d\n', K);
    x_star = ones(n,1);
    x_star(n-K+1:2:end) = -1;
    for m=1:n
        A = normrnd(0,1,m,n,runNumber); 
        y = zeros(m,runNumber);
        xhat = zeros(n,runNumber);
        label = zeros(runNumber,1);
        parfor run = 1:runNumber
            %fprintf('Running %d\n', run);
            y(:,run) = A(:,:,run)*x_star;
            xhat(:,run) = constrainedTV(A(:,:,run), y(:,run));
            if sum(abs(xhat(:,run)-x_star))<thresh
                label(run) = 1;
            end
        end
        prob(m,K+1) = sum(label);
        if prob(m,K+1)==runNumber
            prob(m+1:end,K+1) = runNumber;
            fprintf('Sparisty level %d sample %d, reach full recovery\n', K,m);
            break;
        end
        fprintf('Sparisty level %d sample %d, probability of recovery is %d%%\n', K,m,prob(m,K+1));
    end
end