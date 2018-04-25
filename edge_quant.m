function [diam start edges] = edge_quant(ver,val)

edges = [ver(val,1) ver(val,2)];
diam = norm(max(edges) - min(edges));
%ver(min(val)+1:max(val)-1,:) = [];
start = min(val);
mid = [ver(start,1) + diam/2 ver(start,2)];