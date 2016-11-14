###############################################################
#The Legendre Polynomials:
#The first n Legendre Polynomials form an orthogonal basis of 
#the space of polynomial functions up to order n -1 
###############################################################

function [L, l] = legendre(n)
  % produce the first n legendre polynomials {1, x, ..} 
  % make use of recursive definition:
  % Lk(x) = (2k-1)/k * x * Lk-1(x) - (k-1)/k * Lk-2(x)
  
  L = zeros(n);
  L(1,n) = 1;
  if n > 1
    L(2,n-1) = 1;
  end
  for k = 3:n
    lp = L(k-1,2:n);
    lp(n) = 0;
    lp2 = L(k-2,:);
    lnew = (2*k-1)/k .* lp - (k-1)/k .* lp2;
    L(k,:) = lnew;
  endfor
  l = L(n,:);
endfunction