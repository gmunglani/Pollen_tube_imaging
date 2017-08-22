function [cstart,cstop] = closest_bound(total,xx,y,startpos,stoppos)

dstart = sqrt(sum(bsxfun(@minus, total, [y(startpos) xx(startpos)]).^2,2));
[cstart cstart] = min(dstart);

dstop = sqrt(sum(bsxfun(@minus, total, [y(stoppos) xx(stoppos)]).^2,2));
[cstop cstop] = min(dstop);

