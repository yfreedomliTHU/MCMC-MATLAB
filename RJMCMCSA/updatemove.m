%Perform Update move
i=t;
%P_AIC=k(i)*(c+1)+c(1+d);     %AIC
newmu=zeros(1,d);
if(k(i)==0)
    mu{i+1}=mu{i};
end
for subnum=1:k(i)     
    for j=1:d
      newmu(1,j)= (min(x(:,j))-step(j)) + ((max(x(:,j))+step(j))-(min(x(:,j))-step(j)))*rand(1,1);
    end
    % UPDATE 
    D1 = D;
    D1(:,d+1+subnum) = RBF(newmu,x);
    P=eye(N)-D*inv(D'*D)*D';
    P1=eye(N)-D1*inv(D1'*D1)*D1';
    rRJSA = 1;
    for j=1:c,
      rRJSA= rRJSA * ((y(:,j)'*P*y(:,j))/(y(:,j)'*P1*y(:,j)))^(N/2); 
    end;
    ARJSA = min(1,rRJSA);
    u=rand();
    if (u<ARJSA)
      mu{i+1}(subnum,:) = newmu;
      D=D1;
    else
      mu{i+1}(subnum,:) = mu{i}(subnum,:);
      D=D;
    end; 
end; 
k(i+1) = k(i); % Don't change dimension.











