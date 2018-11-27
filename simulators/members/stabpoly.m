function spoly = stabpoly(n)
% STABPOLY creates a stable (discrete-time) polynomial of order n,
% i.e. all roots has absolute value less than 1.


% prob of an integrator is 0.10 for the first and 0.01 for all others
nint = (rand(1,1)<0.10)+sum(rand(n-1,1)<0.01);

% prob of repeated roots is 0.05
nrepeated = floor(sum(rand(n-nint,1)<0.05)/2);

% prob of complex roots is 0.5
ncomplex = floor(sum(rand(n-nint-2*nrepeated,1)<0.5)/2);
nreal = n-nint-2*nrepeated-2*ncomplex;

% determine random poles
rep = 2*rand(nrepeated,1)-1;
mag = rand(ncomplex,1);
ang = pi*rand(ncomplex,1);
jay = sqrt(-1);
complex = mag.*exp(jay*ang);
re = real(complex); im = imag(complex);

shrink = .9;
poles = [re+jay*im; re-jay*im;...
         ones(nint,1)*shrink; ...
         rep;rep;2*rand(nreal,1)*shrink-1];

% print the polynomial
spoly = poly(poles);