function q = myQuantile(x,p)
% This is a script to calculate the quantiles of a vector if 
% the stats toolbox is not present
% It has been checked against Matlab's own function

N = numel(x);
xx = reshape(sort(x),N,1);

pp = reshape(((1:N)-0.5)/N,N,1);

szout = size(p);

p = reshape(p,numel(p),1);

q = interp1q(pp,xx,p); 

q = reshape(q,szout);


