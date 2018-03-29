%Perform Split move
i=t;
zetastar = 0.1;                            % Merge&split parameter.
Pos = unidrnd(k(i),1,1);     % Choose a centre randomly.
mu0 = mu{i}(Pos,:);
ums = rand(size(mu0));
mu1 = mu0 - ums*zetastar;
mu2 = mu0 + ums*zetastar;
% UPDATE D and P:
tempPos= d+1+Pos;
if (tempPos==d+1+k(i))
  tempD = D(:,1:tempPos-1);    
else
  tempD = [D(:,1:tempPos-1) D(:,tempPos+1:k(i)+d+1)];      
end
D1 = [tempD RBF(mu1,x) RBF(mu2,x)];
P=eye(N)-D*inv(D'*D)*D';
P1=eye(N)-D1*inv(D1'*D1)*D1';
rsplit= zetastar*inv(k(i)+1) * k(i) * exp(-C) * ((y(:,1)'*P*y(:,1))/(y(:,1)'*P1*y(:,1)))^(N/2);      
for j=2:c
  rsplit= rsplit * ((y(:,j)'*P*y(:,j))/(y(:,j)'*P1*y(:,j)))^(N/2); 
end
Asplit=min(1,rsplit);
distance1 = zeros(k(i),1);                 % Calculate Euclidean distance.
distance2 = norm(mu1-mu2); 
for j=1:k(i)
  if (j== Pos)
    distance1(j) = inf; 
  else
    distance1(j)=norm(mu1-mu{i}(j,:));     % Euclidean distance;
  end
  if(distance1(j) < distance2)             % Dont accept;
    Asplit=0;
  end
end; 
u=rand();
if(u<=Asplit)
    tempmu=mu{i};
    if (Pos==1)
      muT = tempmu(2:k(i),:); 
    elseif (Pos==k(i))
      muT = [tempmu(1:k(i)-1,:)];
    else
      muT = [tempmu(1:Pos-1,:); tempmu(Pos+1:k(i),:)];
    end
    mu{i+1} = [muT; mu1; mu2];
    k(i+1)=k(i)+1;
    D=D1;
else
    mu{i+1} = mu{i};
    k(i+1) = k(i);
    D=D;
end


