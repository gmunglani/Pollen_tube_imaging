function [diam start edges] = edge_quant(ver,val)

edges_arr = [ver(val,1) ver(val,2)];
[tmp edge_max] = max(edges_arr(:,1));
[tmp edge_min] = min(edges_arr(:,1));
edges = [edges_arr(edge_max,:); edges_arr(edge_min,:)];
diam = norm(max(edges) - min(edges));
%ver(min(val)+1:max(val)-1,:) = [];
start = min(val);
mid = [ver(start,1) + diam/2 ver(start,2)];