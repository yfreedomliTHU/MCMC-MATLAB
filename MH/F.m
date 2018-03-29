function result = F(x,sigma )
     result=1/(2*pi*sqrt(det(sigma)))*exp(-1/2*(x'/sigma*x));
end
