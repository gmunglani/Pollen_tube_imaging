function [cstart,cstop] = closest_bound(total,xctf,yctf,startpos,stoppos)

xc = [xctf(startpos) xctf(stoppos)];
yc = [yctf(startpos) yctf(stoppos)];

% Calculate the gradient of the center line to get the normals
dx = gradient(xc); dx(find(dx == 0)) = 0.01;
dy = gradient(yc); dy(find(dy == 0)) = -0.01;

nfitc1 = fit(vertcat(xc(1),(xc(1) - dy(1))),vertcat(yc(1),(yc(1) + dx(1))),'poly1');
edge1 = total(:,1) - nfitc1.p1.*total(:,2) - nfitc1.p2;
[tmp cstart] = min(abs(edge1));

nfitc2 = fit(vertcat(xc(2),(xc(2) - dy(2))),vertcat(yc(2),(yc(2) + dx(2))),'poly1');
edge2 = total(:,1) - nfitc2.p1.*total(:,2) - nfitc2.p2;
[tmp cstop] = min(abs(edge2));
