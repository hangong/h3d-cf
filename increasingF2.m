function fout = increasingF2(in,out,p,a)

%% in: vector of x values

%% out: target vector of y values

%% p: small fractional number (e.g. 0.000001)

%% r=z where z=f(x) where f is monotnoically increasing and z~=y

if ~exist('p','var'), p = 0; end
if ~exist('a','var'), a = 255; end
 
%quantise (necessary as we solve per quantisation level)
in = round(in.*a)./a;
 
%We solve for f() by qunatising and histogramming
cross = zeros(a+1,1);
ticks = linspace(0,1,a+1);
edges = linspace(0-1/a*0.5,1+1/a*0.5,a+2);
[w,~,idx] = histcounts(in,edges);

cvv = diag(ticks.^2.*w);
for i = 1:length(out)
    cross(idx(i)) = cross(idx(i)) + out(i);
end
cross = -(ticks').*cross;
 
%We use a difference function to make sure the function is increasing
D_n = zeros(a,a+1);
for i = 1:a
    D_n(i,i:i+1)=[-1 1];
end
D_n_1 = D_n(1:a-1,1:a);
 
%Now we make the 2nd derivative (note the quantised values are
%not evenly split across the domain of x. So we take this into account
%
T = D_n_1*D_n*diag(ticks)*(a.^2);
 
% mult are the per quantisation level scalar such that z=mult*x;
sm = T'*T;
options = optimset('Display','off');
mult = quadprog(cvv+sm.*p,cross,...
    [-D_n*diag(ticks);diag(ticks)],...
    [zeros(a,1);ones(a+1,1)],[],[],zeros(a+1,1),[],[],options);


fin = linspace(0,1,a+1);
fout = fin'.*mult;
