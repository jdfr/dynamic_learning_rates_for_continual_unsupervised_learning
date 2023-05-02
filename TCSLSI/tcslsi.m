function [U, X] = tcslsi(data, varargin)
%TCSLSI Sequrntial algorithm of dimension reduction
%
% [U, X] = TCSLSI(bsize, quality) Uses a sequential algoriritm 
% for dimension reduction of matrix data in stationary case (not dYnamic)
% 
% Input and Output arguments
% data - matrix whose colums are objects and row are attributes
% bsize - size of block
% quality - approximation quality
%
% X - matrix of reduced dimensionality
% U - matrix of factors 
%
% See also TCRDIM
%
% Contributed to TCLIB 
% Copyright (c) by Vitaly Vasilev


%
% Checking Data Parameters
%

A = data;
[m, ndata] = size(A);
Anorm = 0;

U = []; 
X = [];

pnames = {'quality', 'bsize', 'isfullpr', 'iscompr'};
dflts =  {0.25, 10000, 1, 0};
[eid, errmsg, quality, bsize, isfullpr, iscompr] = tcgetargs(pnames, dflts, varargin{:});
if ~isempty(eid)   
    error(sprintf('В tcslsi в процедуре tcdoc_lexicon_rml: %s',eid), errmsg);
end

%
% FIRST LOOP on data
%

% Get first block of data B of size bsize
[B, curpos] = TCGetBlock(A, bsize, 1);
if isempty(B)
    return;
end

Bnorm = sum(sum(B.*B));
Anorm = Anorm + Bnorm; 

% Find eigen decomposition decomposotion of first block
[BV,BS] = tclaneig(@(x) TCBtBfunc(B, x), size(B,2), quality, Bnorm,[]);
BV = BV*diag(sqrt(diag(BS)).^-1);

% Calculate U =  B*BV and make it sparse
U = TCSparseMul(B,BV);

% Check size of data
if curpos == -1
    X = full(BS*BV');
    return;
elseif (isfullpr == 0)
    X = full(BS*BV');
end

while curpos ~= -1    
        
    % Get Next block of data
    [B, curpos] = TCGetBlock(A, bsize, curpos);
    
    if isempty(B)
        break;
    end
       
    % Updating Anorm
    Bnorm = sum(sum(B.*B));
    Anorm = Anorm + Bnorm; 
    
    % Calculate projection of current block to the factor space BP
    BP = U' * B;
    
    % Update matrix X
    if (isfullpr == 0) 
        X = [X, BP];
    end
    
    % Calculate approximation quality of B
    Bquality = sum(sum(BP.*BP)) / Bnorm;
    
    % Calculate desired approximation quality of C = (I-U*U'*B)
    Cdesired = quality - Bquality;
    
    % Check quality level
    if Cdesired > 0.1*quality

        % Find SVD decomposition of matrix C
        [BV,BS] = tclaneig(@(x)TCCtCfunc(B, BP, x), size(B,2),Cdesired,Bnorm,[]);      
        
        % Calculate matrix of factors
        BV = BV*diag(sqrt(diag(BS)).^-1);
                
        % Calculate UB =  B*BV - U*(BP*BV) and make it sparse
        UB = sparse([], [], [], size(B,1), size(BV,2), 2500);
               
        for i=1:size(BV,2)
            ub = B*BV(:,i) - U*(BP*BV(:,i));
            
            [tempu, perm] = sort(ub.*ub);
            tempu = cumsum(tempu / sum(tempu));
            ub(perm(tempu < 0.005)) = 0;             
            
            UB(:,i) = ub;            
        end
                
        % Upate matrix U              
        U = [U, UB]; UB = [];
        
        % Update matrix X
        if (isfullpr == 0) 
            X = [X; zeros(size(BV,2),size(X,2)-size(B,2)), BS*BV'];
        end
                
    end
    
    B = []; BP = [];
end

%
% Дополнительное сжатие данных
%

if iscompr

    % Calculate full projection of data matrix Y = X*X', where  X = U' * A;
    if (isfullpr == 1)
        [B, curpos] = TCGetBlock(A, bsize, 1);
        BP = U'*B;
        Y = BP*BP';
        while curpos ~= -1
            [B, curpos] = TCGetBlock(A, bsize, curpos);
            if isempty(B)
                break;
            end
            BP = U'*B;
            Y = Y + BP*BP';
        end
        Y = full(Y);
    else
        Y = full(X*X');
    end

    % Calculating eigendecomposition of X'*X
    Ynorm = sum(sum(Y));
    [Q, D] = eig(Y);
    [Dsorted, perm] = sort(-diag(D));

    % Selecting eigenvectors
    sel = (cumsum(-Dsorted) / Anorm) <= quality;
    sel(min(sum(sel)+1, length(sel))) = logical(1);
    perm = perm(sel);
    Q = Q(:,perm);

    % Make U = U*Q(:,perm) sparse and minimize memory usage
    U = TCSparseMul(U,Q);

end

%
% Calculating resulting matrix X = Q'*X;
%
if nargout > 1
    X = full(U'*A);
end

%% Локальные функции

function y=TCBtBfunc(B, x)
% TCBTBFUNC Calculating product of B'*B*x
% y=TCBtBfunc(x)
% Function defining a linear operator applied to x. 
%
% Copyright Vasilev Vitaly (2002) (0.45 s)
%

y = ((B*x)'*B)'; % NEW


function y=TCCtCfunc(B, BP, x)
% TCCTCFUNC Calculating product of C'*C*x
% y=TCCtCfunc(x)
% Function defining a linear operator C'*C applied to x, where
% C defined as C = (I-U*U')*B  --> C'*C*x = B'*(B*x) - BP'*(BP*x);
%
% Copyright Vasilev Vitaly (2002) (0.82 s)
%
y = ((B*x)'*B - (BP*x)'*BP)'; % NEW

function [B, curpos] = TCGetBlock(A, bsize, curpos)
% TCGETBLOCK Get Next block of data
% [B, curpos] = TCGetBlock(bsize, curpos) Function load next
% block of data from the storage and return it in matrix B
if min(size(bsize)) > 1
    B = [];
    block_set = [];    
    while length(block_set) == 0
        if curpos == -1
            return;
        end
        block_set = find(bsize(curpos, :));
        curpos = curpos + 1;
        if curpos > size(bsize, 1)
            curpos = -1;
        end
    end
    B = A(:, block_set);
else
    % Обработка блоков фиксированной длины
    B = A(:, curpos:min([(curpos+bsize-1), size(A,2)]) );
    curpos = curpos + bsize;
    if curpos > size(A, 2)
        curpos = -1;
    end
end

function U = TCSparseMulNew(A,B)
% TCSparseMul Calculating U = A*B and make U sparse
U = sparse([], [], [], size(A,1), size(B,2), 5000);
cut = sqrt(0.5/size(A,1));
for i=1:size(B,2)
    ucol = A*B(:,i);
    ucol(abs(ucol) < cut) = 0;
    U(:,i) = ucol;
end

function U = TCSparseMul(A,B)
% TCSparseMul Calculating U = A*B and make U sparse
U = sparse([], [], [], size(A,1), size(B,2), 5000);
for i=1:size(B,2)
    ucol = A*B(:,i);
    [tempu, permu] = sort(ucol.*ucol);
    tempu = cumsum(tempu / sum(tempu));   
    ucol(permu(tempu < 0.001)) = 0;        
    U(:,i) = ucol;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequential variant of X = U' * A
% % Get first block of data B of size bsize
% [B, curpos] = TCGetBlock(bsize, 1);
% 
% % Calculate projection of block on full space
% BP = U'*B;
% 
% % Calculate XX
% XX  = BP*BP';
% X = BP;
% while curpos < ndata
%     
%     % Get Next block of data
%     [B, curpos] = TCGetBlock(bsize, curpos);    
%     
%     % Calculate projection of block on full space
%     BP = U'*B;
%     
%     % Updating X
%     X = [X, BP]; 
%     
%     % Upadate XX
%     XX = XX + BP*BP';    
% end