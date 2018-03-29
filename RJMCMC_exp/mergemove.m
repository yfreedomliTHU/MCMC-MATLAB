%Perform Merge move
i=t;
zetastar = 0.1;                            % Merge&split parameter.
Pos = unidrnd(k(i),1,1);     % Choose a centre randomly.
mu1 = mu{i}(Pos,:);
distance = zeros(k(i),1);
for j=1:k(i)
  if (j== Pos)
    distance(j) = inf; 
  else
    distance(j)=norm(mu1-mu{i}(j,:));      % Euclidean distance;
  end
end
Pos2 = find(distance == min(distance));
mu2 = mu{i}(Pos2,:);
mu0 = 0.5*(mu1 + mu2);
%Update D&mu
%Delete mu1
tempPos1 = 1+d+Pos;
if (tempPos1==d+1+k(i))
  tempD = D(:,1:tempPos1-1);     
else
  tempD = [D(:,1:tempPos1-1) D(:,tempPos1+1:k(i)+d+1)];      
end
%Delete mu2
tempPos2 = Pos2+d+1;
if (Pos2>Pos)  %Attention:mu1 has been deleted. 
  if (tempPos2==d+1+k(i))
    tempD = tempD(:,1:tempPos2-2);     
  else
    tempD = [tempD(:,1:tempPos2-2) tempD(:,tempPos2:k(i)+d)];      
  end
elseif (Pos2<Pos)
    tempD = [tempD(:,1:tempPos2-1) tempD(:,tempPos2+1:k(i)+d)];      
else
    error('Merge move wrong!');
end
D1 = [tempD RBF(mu0,x)];
P=eye(N)-D*inv(D'*D)*D';
P1=eye(N)-D1*inv(D1'*D1)*D1';
rmerge= (k(i) * exp(C))/(zetastar*(k(i)-1))* ((y(:,1)'*P*y(:,1))/(y(:,1)'*P1*y(:,1)))^(N/2);      
for j=2:c
  rmerge= rmerge * ((y(:,j)'*P*y(:,j))/(y(:,j)'*P1*y(:,j)))^(N/2); 
end
Amerge=min(1,rmerge);
if (min(distance)<2*zetastar)   % To ensure reversibility.
   Amerge = 0;   
end;
u=rand();
if(u<=Amerge)
    tempmu=mu{i};
    if (Pos==1)
      muT = tempmu(2:k(i),:); 
    elseif (Pos==k(i))
      muT = [tempmu(1:k(i)-1,:)];
    else
      muT = [tempmu(1:Pos-1,:); tempmu(Pos+1:k(i),:)];
    end
    if (Pos2>Pos)
       if (Pos==k(i))
           muT = muT(1:k(i)-2,:);
       else
           muT = [muT(1:Pos2-2,:); muT(Pos2:k(i)-1,:)];
       end
    elseif (Pos2<Pos)
       if (Pos2==1)
           muT = muT(2:k(i)-1,:); 
       else
           muT = [muT(1:Pos2-1,:); muT(Pos2+1:k(i)-1,:)];
       end
    else
    error('Merge move wrong!');
    end
    mu{i+1} = [muT; mu0];
    k(i+1) = k(i)-1;
    D=D1;
else
    mu{i+1} = mu{i};
    k(i+1) = k(i);
    D=D;
end




