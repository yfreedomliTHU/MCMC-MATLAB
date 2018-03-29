clear all,close all,clc,tic;  %MCMCRBF
load data1.mat;
%load data2.mat;
[N,d]=size(x);             % Get dimension and num of x.
c=2;                       % Dimension of y.
Simulation_time = 1000;    % Length of the Markov chain simulation.
k = ones(Simulation_time,1)*20;   % Model order(number of basis).
mu = cell(Simulation_time,1); % Radial basis centres.
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
for i=1:Simulation_time-1
   %P_AIC=k(i)*(c+1)+c(1+d);     %AIC
   newmu=zeros(1,d);
   for subnum=1:k(i)     
      for j=1:d
         newmu(1,j)= (min(x(:,j))-step(j)) + ((max(x(:,j))+step(j))-(min(x(:,j))-step(j)))*rand(1,1);
      end
      D1 = D;
      D1(:,d+1+subnum) = RBF(newmu,x);
      P=eye(N)-D*inv(D'*D)*D';
      P1=eye(N)-D1*inv(D1'*D1)*D1';
      r = 1;
      for j=1:c,
         r= r * ((y(:,j)'*P*y(:,j))/(y(:,j)'*P1*y(:,j)))^(N/2); 
      end;
      AR = min(1,r);
      u=rand();
      if (u<AR)
         mu{i+1}(subnum,:) = newmu;
         D=D1;
      else
         mu{i+1}(subnum,:) = mu{i}(subnum,:);
      end; 
   end 
    k(i+1) = k(i); % Don't change dimension.
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
v1=Dtest*alpha;
save('MCMCRBF2015013178.mat','v1');
%save('MCMCRBF2015013178.mat','v2','-append');
toc;