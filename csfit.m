function pp = csfit(x,y,n)

breaks = linspace(0,1,n+1);
[xcount,xind] = histc(x,[0,breaks(2:end-1),1]);

w = xcount/sum(xcount);

pointMean = zeros(n,2);
for i = 1:n
    pointMean(i,1) = nanmean(x(xind==i));
    pointMean(i,2) = nanmean(y(xind==i));
end

pointMean = pointMean(w>0,:);

pp = pchip(pointMean(:,1),pointMean(:,2));

%fnplt(pp); hold on;
%plot(pointMean(:,1),pointMean(:,2),'+');
%error('h');
