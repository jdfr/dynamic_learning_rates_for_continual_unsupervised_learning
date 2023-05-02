function [V,D, bnd,j,work, res_quality] = TCLANEIG(A,n,quality, Anorm, options)

%TCLANEIG  Compute factor space for specified approx. quality.
%   TCLANEIG uses lanczos algorithm from PROPACK for computing
%   optimal linear factor space that guranty specified approximation 
%   quality of real simmetric matrix A. 
%
%   Calling sequence is
%
%   [V,D,ERR] = LANEIG(A,QUALITY,ANORM,OPTIONS)
%   [V,D,ERR] = LANEIG(@Afun,N,QUALITY,ANORM,OPTIONS)
%
%   Afun - function that realize linear operator A*x
%   n      - dimension of square matrix A
%   quality - relative quality of approximation of A in frobenius norm
%   Anorm - square of frobenius norm of matrix A
%   
%   On exit ERR contains the computed error bounds.  K is the number of
%   eigenvalues desired and SIGMA is numerical shift or a two letter string
%   which specifies which part of the spectrum should be computed:
%
%   The OPTIONS structure specifies certain parameters in the algorithm.
%
%    Field name      Parameter                              Default
%   
%    OPTIONS.tol     Convergence tolerance                  16*eps
%    OPTIONS.lanmax  Dimension of the Lanczos basis.
%    OPTIONS.v0      Starting vector for the Lanczos        rand(n,1)-0.5
%                    iteration.
%    OPTIONS.delta   Level of orthogonality among the       sqrt(eps/K)
%                    Lanczos vectors.
%    OPTIONS.eta     Level of orthogonality after           10*eps^(3/4)
%                    reorthogonalization. 
%    OPTIONS.cgs     reorthogonalization method used        0
%                    '0' : iterated modified Gram-Schmidt 
%                    '1' : iterated classical Gram-Schmidt
%    OPTIONS.elr     If equal to 1 then extended local      1
%                    reorthogonalization is enforced. 
%    
%
%   See also LANPRO, EIGS, EIG.

% References: 
% R.M. Larsen, Ph.D. Thesis, Aarhus University, 1998.
%
% B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
% Prentice-Hall, Englewood Cliffs, NJ, 1980.
%
% H. D. Simon, ``The Lanczos algorithm with partial reorthogonalization'',
% Math. Comp. 42 (1984), no. 165, 115--142.

% This file is corrected to using function handles by Vasiljev Vitaly (2002)
% for finding eigen values that guranty specified  approximation quality of matrix A.

% Rasmus Munk Larsen, DAIMI, 1998


%%%%%%%%%%%%%%%%%%%%% Parse and check input arguments. %%%%%%%%%%%%%%%%%%%%%%
if nargin<3
    error('TCLANEIG', 'В TCLANEIG недостаточно входных параметров');
end

% Quick return for n<2  or quality <= 0
% set k =1
if n < 1 |  quality <= 0
    if nargout < 2
        V = zeros(1,1);
    else
        V = eye(n,1);
        D = zeros(1,1);
        bnd =zeros(1,1);
    end
    return
end

if n == 1 
    D = feval(A,1);
    V = 1;
    dnb = 0;
    if nargout<2
        V=D;
    end
    return
end

% Estimating starting parameters: k, lanmax, end_k, end_lanmax
k = min([15, n]);
lanmax	=	2*k + 2;
tol = 16*eps;  % here were 16*eps
r = rand(n,1)-0.5;
end_lanmax = min(600,n);

% Parse options struct for user parameters
if ~isempty(options) & isstruct(options)
    c = fieldnames(options);
    for i=1:length(c)
        if strmatch(c(i),'v0'), r = getfield(options,'v0'); r=r(:); end
        if strmatch(c(i),'tol'), tol = getfield(options,'tol'); end
        if strmatch(c(i),'lanmax'), lanmax = getfield(options,'lanmax'); end        
        if strmatch(c(i),'end_lanmax'), end_lanmax = getfield(options,'end_lanmax'); end        
    end
end

% Protect against absurd arguments.
tol = max(tol,eps);
lanmax = min(lanmax,n);
end_lanmax = min(end_lanmax, n); 
if size(r,1)~=n
    error('MATLAB:TCLANEIG', 'В TCLANEIG v0 должен быть вектором длины n')
end

if k>lanmax
    error('MATLAB:TCLANEIG', 'В TCLANEIG должно K <= LANMAX <= N.');
end
ksave = k;

neig = 0; nrestart=-1;
j = min(2*k+2,lanmax);

%%%%%%%%%%%%%%%%%%%%% Here begins the computation  %%%%%%%%%%%%%%%%%%%%%%

V = []; T = []; anorm = []; work = zeros(1,2); rnorm=-1;

while neig < k 
    
    %%%%%%%%%%%%%%%%%%%%% Compute Lanczos tridiagonalization %%%%%%%%%%%%%%%%%
    j = min(lanmax,j+1-mod(j,2));
    
    % "Trick" to avoid unwanted zero eigenvalues when laneig is used for
    % SVD calculations. (Nothing to if lanmax is odd, though.)
    % j - contain number of lanczos steps
    [V,T,r,anorm,ierr,w] = lanpro(A,n,j,r,options,V,T,anorm);
    
    work= work + w;
    
    if ierr<0 % Invariant subspace of dimension -ierr found. 
        j = -ierr;
    end
    
    %%%%%%%%%%%%%%%%%% Compute eigenvalues and error bounds %%%%%%%%%%%%%%%%%%
    
    % Analyze T
    [D,top,bot,err] = tqlb([full(diag(T))],full([0;diag(T,1)]));
    
    [D,I] = sort(D);
    bot = bot(I);
    
    % Set simple error bounds
    rnorm = norm(r);
    bnd = rnorm*abs(bot);
    
    % Use Largest Ritz value to estimate ||A||_2. This might save some
    % reorth. in case of restart.
    anorm = max(abs(D));
    
    % Estimate gap structure and refine error bounds
    bnd = refinebounds(D,bnd,n*eps*anorm);
    
    %%%%%%%%%%%%%%%%%%% Check convergence criterion %%%%%%%%%%%%%%%%%%%%
    % Reorder eigenvalues according to SIGMA
    
    [dummy,IPART] = sort(-abs(D));
    
    D = D(IPART);  bnd = bnd(IPART);
    
    % Check if enough have converged.
    neig = 0; % contain number of converged singular values
    for i=1:min(j,k)
        if bnd(i) <= tol*abs(D(i))
            neig = neig + 1;
        end
    end
    if neig < 1
        neig = 1;
    end
    
    %%%%%%%%%%% Check whether to stop or to extend the Krylov basis? %%%%%%%%%%
              
    if ierr<0 % Invariant subspace found
        if j<k
            warning(['Invariant subspace of dimension ',num2str(j-1),' found.'])
        end
        break;    % terminate loop and take j as k
    end

    cur_quality = cumsum(abs(D(1:neig))) / Anorm;
    
    if (max(cur_quality) >= quality)  |  (j >= end_lanmax)
        
        % Select minimun number of k
        k = sum(cur_quality < quality) + 1;
        
        % Save result quality
        try
            res_quality = cur_quality(min(k,neig));
        catch
            disp('laneig catch');
        end

        break;  % terminate loop        
        
    elseif neig < k % Number of converged eigen values is small       
    
        j = j + 1 + max(10, 0.5*(k-neig)*j/(neig+1));
        if j >= lanmax 
            lanmax = end_lanmax;
            j = min(lanmax,j);
        end        
        
    else % Estimate new value for k using regression model y = a*(x^b)
    
        % Find regression coefficients using last 10 eigenvalues
        seqint = (neig-min(neig, 10)+1):neig;
        cumseq = cumsum(abs(D(1:neig)))/Anorm;       
        [a, b] = regress_axb(cumseq(seqint), seqint');        
        
        % Estimate appropriate new value for k
        new_k = neig+1; 
        while (a*new_k^b < quality) & new_k < end_lanmax
             new_k = new_k +1;
        end
        
        k = new_k;              
        j = j + 1 + max(15, 0.5*(k-neig)*j/(neig+1));               
        if j >= lanmax 
            lanmax = min([end_lanmax, j*2+2]);
            j = min(lanmax,j);
        end                
    end     
   
    
    nrestart = nrestart + 1;    
end


%%%%%%%%%%%%%%%% Lanczos converged (or failed). Prepare output %%%%%%%%%%%%%%%
k = min(k,j);

if nargout>1
    j = size(T,1);
    [Q,D] = eig(full(T)); D = diag(D);
    [D,I] = sort(D);
    % Compute and normalize Ritz vectors (overwrite V to save memory).
    V = V*Q(:,I(IPART(1:k)));
    for i=1:k
        nq = norm(V(:,i));
        if isfinite(nq) & nq~=0 & nq~=1
            V(:,i) = V(:,i)/nq;
        end
    end
    [D,I] = sort(D);
    D = D(IPART(1:k));
end

% Pick out desired part of the spectrum
if length(D)~=k
    D = D(1:k);
    bnd = bnd(1:k);
end

if nargout<2
    V = D;
else
    D = diag(D);
end


function [a, b] = regress_axb(y, x)
lx = log(x);
ly = log(y);
my = mean(ly);
mx = mean(lx);
myx = (ly'*lx)/length(y);
mxx = (lx'*lx)/length(x);
b = (myx - my * mx) / (mxx - mx*mx);
a = exp(my - mx*(myx-my*mx)/(mxx-mx*mx));