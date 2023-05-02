function [lambda,top,bot,err] = tqlb_matlab(alpha,beta)

% TQLB: Compute eigenvalues and top and bottom elements of
%       eigenvectors of a symmetric tridiagonal matrix T.
%
% [lambda,top,bot,err] = tqlb_matlab(alpha,beta)
%
% Input parameters:
%   alpha(1:n)   : Diagonal elements.
%   beta(2:n)    : Off-diagonal elements.
% Output parameters:
%   lambda(1:n)  : Computed eigenvalues.
%   top(1:n)     : Top elements in eigenvectors.
%   bot(1:n)     : Bottom elements in eigenvectors.
%   err          : dummy argument.


% Rasmus Munk Larsen, DAIMI, 1998


%
% This is a slow Matlab substitute for the 
% TQLB MEX-file.
%

n = length(alpha);
T = spdiags([[beta(2:n);0] alpha(1:n) beta(1:n)],-1:1,n,n);

safe_and_slow = 1;
if safe_and_slow==1
  [V,lambda] = eig(full(T)); lambda = diag(lambda);
  bot = V(end,:)';
  top = V(1,:)';
  err=0;
else
  lambda = eig(T);
  bot  = paige(T,lambda); 
  % Use Paiges formula for bottom elements. May result in overflow 
  % and inaccurate error bounds for large n;
  top = [];
  err=0;
end

