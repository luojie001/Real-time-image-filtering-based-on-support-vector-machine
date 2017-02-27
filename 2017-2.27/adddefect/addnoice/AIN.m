clc;clear;
lambda=5;
sw=3;
sl=3;
sizeOut=[sw,sl];
numelOut=sw*sl;
lambda=ones(numelOut,1)*lambda;
r=zeros(sw,sl);

% For large lambda, use the method of Ahrens and Dieter as
% described in Knuth, Volume 2, 1998 edition.
k = find(15 <= lambda & lambda < Inf);
if ~isempty(k)
   alpha = 7/8;
   lk = lambda(k);
   m = floor(alpha * lk);

   % Generate m waiting times, all at once
   x = randg(m);
   t = (x <= lk);

   % If we did not overshoot, then the number of additional times
   % has a Poisson distribution with a smaller mean.
   r(k(t)) = m(t) + poissrnd(lk(t)-x(t));

   % If we did overshoot, then the times up to m-1 are uniformly
   % distributed on the interval to x, so the count of times less
   % than lambda has a binomial distribution.
   if ~all(t)
       r(k(~t)) = binornd(m(~t)-1, lk(~t)./x(~t));
   end
end
% For small lambda, generate and count waiting times.
j = find(lambda < 15);
p = zeros(numel(j),1,'like',lambda);
while ~isempty(j)
    p = p - log(rand(numel(j),1,'like',lambda));
    t = (p < lambda(j));
    j = j(t);
    p = p(t);
    r(j) = r(j) + 1;
end
r