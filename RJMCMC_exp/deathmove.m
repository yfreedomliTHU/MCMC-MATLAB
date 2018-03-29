% Perform Death move
i=t;
deletemuPos= d+1+unidrnd(k(i),1,1);    %Choose the basis centre to be deleted
if (deletemuPos==d+1+k(i))
  D1 = [D(:,1:deletemuPos-1)];     
else
  D1 = [D(:,1:deletemuPos-1) D(:,deletemuPos+1:k(i)+d+1)];      
end
P=eye(N)-D*inv(D'*D)*D';
P1=eye(N)-D1*inv(D1'*D1)*D1';
rdeath= k(i) * exp(C) * ((y(:,1)'*P*y(:,1))/(y(:,1)'*P1*y(:,1)))^(N/2);      
for j=2:c
  rdeath= rdeath * ((y(:,j)'*P*y(:,j))/(y(:,j)'*P1*y(:,j)))^(N/2); 
end
u=rand();
Adeath=min(1,rdeath);
if(u<=Adeath)
    tempmu=mu{i};
    if(deletemuPos==1+d+1)
         mu{i+1} = tempmu(2:k(i),:); 
    elseif(deletemuPos==(1+d+k(i)))
         mu{i+1} = tempmu(1:k(i)-1,:);
    else
         mu{i+1} = [tempmu(1:deletemuPos-1-d-1,:); tempmu(deletemuPos-d-1+1:k(i),:)];
    end
    k(i+1) = k(i)-1;
    D=D1;
else
    mu{i+1} = mu{i};
    k(i+1) = k(i);
    D=D;
end