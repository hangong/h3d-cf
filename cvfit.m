function f = cvfit(x,y,method)

if ~exist('method','var'), method = 'quad'; end

msk = x>0 & y>0;

switch method
    case 'quad'
        f = increasingF2(x(msk),y(msk),1e-5);
    case 'hist'
        h = imhist(y(msk));
        [~,f] = histeq(x(msk),h);
        f = f';
end

% plot(x*1000,y,'.');  hold on;
% plot([0:1000]',f);
% error('h');
