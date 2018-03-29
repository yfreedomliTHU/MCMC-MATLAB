clear all,close all,clc,tic;  % Uniform distribution
mu=[5;10];sigma=[1,1;1,4];
simulation_times=5*1e5;       % Initial simulation_time
get_num=min(simulation_times/50,1000);
S_size=simulation_times;
S=zeros(2,S_size);
S(1,:)=unifrnd(-1,1,1,simulation_times)*5;
S(2,:)=unifrnd(-1,1,1,simulation_times)*10;
X=[];
Y1=F([0;0],sigma);
for i=1:simulation_times
    Y=F(S(:,i),sigma);
    u=rand();
    r=min(1,Y/Y1);
    if r>u
        Y1=Y;
        X=[X S(:,i)+mu];
    end  
end
X1=X(:,get_num:end);
% plot
theorical_data=mvnrnd(mu,sigma,simulation_times-get_num);
subplot(1,2,1),plot(theorical_data(:,1),theorical_data(:,2),'.'),
set(gca,'XLim',[0,10],'YLim',[0,20]),
title('理论取值','FontSize',14);
subplot(1,2,2),
plot(X1(1,:),X1(2,:),'.'),
set(gca,'XLim',[0,10],'YLim',[0,20]),
title('MH采样结果','FontSize',14);
figure;subplot(1,2,1);
hist3(theorical_data,[100,200]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('理论取值');
subplot(1,2,2);
hist3(X1',[100,200]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('MH采样结果');
% calculate the correlation coefficient
corrmatrix=corrcoef(X1(1,:),X1(2,:));
p=corrmatrix(2,1)
e=vpa(100*(abs(p-0.5)/0.5),4);  % calculate error
disp(['e=',char(e),'%'])
toc;