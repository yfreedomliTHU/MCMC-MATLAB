function [ y ] = Gaussian_noise( sigma,x)
     [a,b]=size(x);
     [c,d]=size(sigma);
     y = zeros(a,c);
     for i=1:a
         for j=1:c
             y(i,j)= 1/sqrt(2*pi*sigma(j))*exp(-1*x(i,j)^2/(2*sigma(j)));
         end
     end
end

