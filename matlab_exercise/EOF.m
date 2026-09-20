function [L, EOFs, EC, error] = EOF( U, n, varargin )
% EOF - computes EOF of a matrix.
%
% Usage: [L, EOFs, EC, error] = EOF( M, num, ... )
%
% M is the matrix on which to perform the EOF.  num is the number of EOFs to
% return.  If num='all', then all EOFs are returned.  This is the default.
%
% ... are extra arguments to be given to the svds function.  These will
% be ignored in the case that all EOFs are to be returned, in which case
% the svd function is used instead. Use these with care.
%
% Data is not detrended before handling.  Use the detrend function to fix
% that. 
%
% L are the eigenvalues of the covariance matrix ( ie. they are normalized
% by 1/(m-1), where m is the number of rows ).  EC are the expansion
% coefficients (PCs in other terminology) and error is the reconstruction
% error (L2-norm).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: EOF.m,v 1.1 2002/10/08 22:23:36 dmk Exp $	
%
% Copyright (C) 2001 David M. Kaplan
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
  n = 'all';
end

s = size(U);
ss = min(s);

if (ischar(n) & n == 'all') | n >= ss
  % Use svd in case we want all EOFs - quicker.
  [ C, lambda, EOFs ] = svd( full(U) ); 
else
  % Otherwise use svds.
  [ C, lambda, EOFs, flag ] = svds( U, n, varargin{:} );
  
  if flag % Case where things did not converge - probably an error.
    warning( 'DMK - Eigenvalues did not seem to converge!!!' );
  end
  
end

EC = C * lambda; % Expansion coefficients.
L = diag( lambda ) .^ 2 / (s(1)-1); % eigenvalues.

% Compute error.
diff=(U-EC*EOFs');
error=sqrt( sum( diff .* conj(diff) ) );

