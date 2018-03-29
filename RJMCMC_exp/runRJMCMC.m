clear all,close all,clc,tic;      %RJMCMC experiment
N=100;c=1;d=1;                    % Initial sample number.
x=10*rand(N,1)-5;  
xtest = -5:0.1:5;
y=1.1*(1-5*x+2*x.^3).*exp(-x.^2/2);
Simulation_time = 500;            % Length of the Markov chain simulation.
k = ones(Simulation_time,1)*20;   % Model order(number of basis).
mu = cell(Simulation_time,1);     % Radial basis centres.
C=c+1;                            % AIC=C+1.
step=0.1*(max(x)-min(x));
for i=1:Simulation_time
    mu{i}=zeros(k(i),d);          % Initial mu and Guaranteed within x
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
    u=rand();
    if(u<=0.2)                     % To simplify the experiment set bk=dk=sk=mk=uk=0.2.
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
ypred=D*alpha;
erro=norm(y-ypred)/norm(y)*100

%TEST================================%
[d1,N1]=size(xtest);
k_test=k(t+1);
Dtest=zeros(N1,1+k_test+d1);
Dtest(:,1)=ones(N1,1);
Dtest(:,2:d1+1) = xtest';
for i=d1+2:k_test+d1+1
     Dtest(:,i) = RBF(mu{t+1}(i-d1-1,:),xtest');
end
v1=Dtest*alpha;   %Predicted Number
%PLOT
figure;
hold on;
grid;
plot(x,y,'k+');
plot(xtest,v1,'b-');
xlabel('x');
ylabel('y');

toc;