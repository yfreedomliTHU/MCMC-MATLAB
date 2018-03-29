%Perform Birth move
i=t;
newmu=zeros(1,d);%Propose a new RBF centre at random
for j=1:d
  newmu(1,j)= (min(x(:,j))-step(j)) + ((max(x(:,j))+step(j))-(min(x(:,j))-step(j)))*rand(1,1);
end
D1=[D RBF(newmu,x)];
P=eye(N)-D*inv(D'*D)*D';
P1=eye(N)-D1*inv(D1'*D1)*D1';
rbirth= inv(k(i)+1) * exp(-C) * ((y(:,1)'*P*y(:,1))/(y(:,1)'*P1*y(:,1)))^(N/2);      
for j=2:c
  rbirth= rbirth * ((y(:,j)'*P*y(:,j))/(y(:,j)'*P1*y(:,j)))^(N/2); 
end
u=rand();
Abirth=min(1,rbirth);
if(u<=Abirth)
    mu{i+1} = [mu{i}; newmu];
    k(i+1) = k(i)+1;
    D=D1;
else
    mu{i+1} = mu{i};
    k(i+1) = k(i);
    D=D;
end