function type = find_orient(M)

maxx = size(M,1); maxy = size(M,2); type = 0;
places = [sum(M(1,:)) sum(M(:,maxy)) sum(M(maxx,:)) sum(M(:,1))];
[temp type] = max(places);