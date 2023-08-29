function [p,w,lam,mn,sd] = pca_arm(x)

%*************************************************************************
%
%   FUNCTION:      pca_arm.m
%   =========      ==========
%
%   DESCRIPTION:   My own function for finding principal components.
%   ============   In addition to using terminology I can understand,
%                  my function normalizes the original data matrix
%                  by subtracting the mean and dividing by the standard
%                  deviation.  
%
%                  Inputs:   x   - m-by-n data matrix
%
%                  Returns structure:  
%                    p   = Transformed data ("Scores")
%                    w   = Transformation matrix ("Coeffs")
%                    lam = Eigenvalues of covariance matrix ("Latent")
%                    mn  = Mean values used in normalization
%                    sd  = Standard deviations used in normalization
%                                      
%                  The original data could be reconstructed from:
%                  x = (p*w')*st + mn on a column-by-column basis                 
%
%   COPYWRITE:     Allan R. Moser
%   ==========     Cira Discovery Sciences Consulting
%                  815 Westdale Ave.
%                  Swarthmore, PA  19081
%
%   DATE CREATED:  25-Apr-2014
%   =============
%
%   LAST CHANGED:  25-Apr-2014
%   =============
%
%**************************************************************************


[m, n] = size(x);
mn = mean(x);
sd = std(x);
for k=1:n
    y(:,k) = (x(:,k) - mn(k))/sd(k);
end

% =========================================================
% This code commented out does pca using the covariance matrix.
% This straight-forward method squares the condition number
% of the matrix, which is a bad thing. See the code below, which does 
% pca on data matrix directly
% c = cov(y);
% [u,s,v] = svd(c);
% lam = diag(s);
% w = v;
% =========================================================
% Better method that doesn't square condition number
[u,s,v] = svd(y);
lam = diag(s).^2/(m-1);
w = v;

% There is an indeterminancy in the signs of the components of columns in
% the transformation matrix. Resolve this by changing signs so that larges
% value in a column has a positive sign
for k=1:n
    [msv,ksv] = max(abs(v(:,k)));
    if (v(ksv,k) < 0)
        w(:,k) = -1*v(:,k);
    end
end
p = y*w;
end
