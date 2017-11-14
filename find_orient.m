function type = find_orient(M)

type = 0;
places = [sum(M(1,:)) sum(M(:,end)) sum(M(end,:)) sum(M(:,1))];
[temp type] = max(places);