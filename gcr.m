% gcr.m

function [x,rho,chi,c,a]=gcr(A,b,n)

b=b(:);

rho(:,1) = b; % first residuum
chi(:,1) = A*rho(:,1);
a(1,1) = 1/sqrt(chi(:,1)'*chi(:,1));
chi(:,1) = chi(:,1)*a(1,1);
c(1) = chi(:,1)'*b;
x = a(1,1)*c(1)*rho(:,1);


for k=2:n,

   % calculate new residuum
   rho(:,k) = rho(:,k-1) - c(k-1)*chi(:,k-1);
   
   chibar = A*rho(:,k);
   
   % orthonormalization
   for j=1:k-1,
      a(k,j) = -chi(:,j)'*chibar;
      chibar = chibar + a(k,j)*chi(:,j);
   end;
   a(k,k) = 1/sqrt(chibar'*chibar);
   chi(:,k) = chibar*a(k,k);
   
   % get the coefficients
   for j=1:k-1,
      a(k,j) = a(k,j)*a(k,k);
   end;
   for j=1:k-1,
      a(k,j) = a(k,j)*a(j,j);
      for l=j+1:k-1,
         a(k,j) = a(k,j) + a(k,l)*a(l,j);
      end;
   end;
   
   % calculate next c
   c(k) = chi(:,k)'*b;
   
   %update solution
   for j=1:k,
      x = x + c(k)*a(k,j)*rho(:,j);
   end;   
end;
