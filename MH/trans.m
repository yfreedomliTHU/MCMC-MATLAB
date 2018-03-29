function result = trans(x,sigma )
     result=1/sqrt(2*pi*det(sigma))*exp(-1/2*(x'/sigma*x));
end
