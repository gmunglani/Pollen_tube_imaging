function [cstart,cstop] = closest_bound(total,xctf,yctf,startpos,stoppos)

xc = [xctf(startpos) xctf(stoppos)];
yc = [yctf(startpos) yctf(stoppos)];

% Calculate the gradient of the center line to get the normals
dx = gradient(xc);
dy = gradient(yc);

nfitc1 = polyfit(vertcat(xc(1),(xc(1) - dy(1))),vertcat(yc(1),(yc(1) + dx(1))),1);
edge1 = total(:,1) - nfitc1(1).*total(:,2) - nfitc1(2);
[tmp cstart] = min(abs(edge1));

nfitc2 = polyfit(vertcat(xc(2),(xc(2) - dy(2))),vertcat(yc(2),(yc(2) + dx(2))),1);
edge2 = total(:,1) - nfitc2(1).*total(:,2) - nfitc2(2);
[tmp cstop] = min(abs(edge2));
        
        
%dstart = sqrt(sum(bsxfun(@minus, total, [y(startpos) xx(startpos)]).^2,2));
%[cstart cstart] = min(dstart);

%dstop = sqrt(sum(bsxfun(@minus, total, [y(stoppos) xx(stoppos)]).^2,2));
%[cstop cstop] = min(dstop);

