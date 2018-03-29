function [result] = RBF(mu,x)    % RBF Function
                                 % Choose Gaussian.
    [n,d] = size(x);             % N = number of data, d = dimension of x.
    result=zeros(n,1);
    for i=1:n
        z=norm(x(i,:)-mu(1,:));  % Euclidean distance.
        %result(i,1)=exp(-(0.2)*z.^(2));%AIC  v1    
        %result(i,1)=5*exp(-(0.2)*z.^(2)); %AIC  v2
        result(i,1)=exp(-(0.3)*z.^(2));%BIC&MDL  v1   
    end
end
