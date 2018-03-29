clear all,close all,clc,tic;%RJSA
%load data1.mat;
load data2.mat;
[N,d]=size(x);             % Get dimension and num of x.
c=2;                       % Dimension of y.
Simulation_time = 100;    % Length of the Markov chain simulation.
T1 = 1;                    % Initial temperature.
T2 = 1e-3;                 % Final temperature.
T_step = (T2-T1)/Simulation_time; % Cooling parameter 
Temp =ones(Simulation_time,1);% SA temperature.
k = ones(Simulation_time,1)*20;   % Model order(number of basis).
mu = cell(Simulation_time,1); % Radial basis centres.
%C=c+1;                        % AIC=C+1.
C=(c+1)*log(N)/2;               % BIC&MDL
step=0.1*(max(x)-min(x));
for i=1:Simulation_time
    mu{i}=zeros(k(i),d);      %Initial mu1 and Guaranteed within x
    for t=1:d             
       mu{i}(:,t)= (min(x(:,t))-step(t))*ones(k(i),1) + ((max(x(:,t))+step(t))-(min(x(:,t))-step(t)))*rand(k(i),1);
    end
end
%Initial D
D=zeros(N,k(1)+d+1);
D(:,1)=ones(N,1);
D(:,2:d+1) = x;
for i=d+2:k(1)+d+1
    D(:,i) = RBF(mu{1}(i-d-1,:),x);
end
rbirth=1;rdeath=1;ratio=1;rmerge=1;rupdate=1;
for t=1:Simulation_time-1
    
    preD=D;
    Temp(t) = T1+T_step*t;          % Cooling schedule. 
    u=rand();
    if(u<=0.2)                  % RJSA set bk=dk=sk=mk=uk=0.2.
        birthmove;
    elseif(u<=0.4)
        deathmove;
    elseif(u<=0.6)
        splitmove;
    elseif(u<=0.8)
        mergemove;
    else
        updatemove;
    end
    %ANNEAL
    P=eye(N)-preD*inv(preD'*preD)*preD';
    P1=eye(N)-D*inv(D'*D)*D';
    Ratio =  exp(-C)*((y(:,1)'*P*y(:,1))/(y(:,1)'*P1*y(:,1)))^(N/2);      
    for j=2:c,
        Ratio = Ratio * ((y(:,j)'*P*y(:,j))/(y(:,j)'*P1*y(:,j)))^(N/2); 
    end;
    rSA=(Ratio)^(1/Temp(t)-1);
    ASA=min(1,rSA);
    u=rand();
    if (u<=ASA)
       mu{t+1} = mu{t+1};
       k(t+1) = k(t+1);
       D=D;
    else
       mu{t+1} = mu{t};
       k(t+1) = k(t);
       D=preD;
    end 
end
sigma=zeros(c,1);
for i=1:c
    sigma(i)=1/N*y(:,i)'*P1*y(:,i);
end
alpha=[];
temp=[];
for j=1:c
    temp=inv(D'*D)*D'*y(:,j);
    alpha=[alpha temp];
end
noise=Gaussian_noise(sigma,x);%Gaussian noise
ypred=D*alpha
erro=norm(y-ypred)/norm(y)*100
%TEST================================%
[N1,d1]=size(xtest);
k_test=k(t+1);
Dtest=zeros(N1,1+k_test+d1);
Dtest(:,1)=ones(N1,1);
Dtest(:,2:d1+1) = xtest;
for i=d1+2:k_test+d1+1
    Dtest(:,i) = RBF(mu{t+1}(i-d1-1,:),xtest);
end
v2=Dtest*alpha;
%save('RJMCMC2015013178.mat','v1');
%save('RJMCMC2015013178.mat','v2','-append');
%save('BIC2015s013178.mat','v1');
%save('BIC2015013178.mat','v2','-append');
toc;