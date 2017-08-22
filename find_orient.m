function type = find_orient(M)

maxx = size(M,1); maxy = size(M,2); type = 0;
if ~isempty(nonzeros(M(1,:))) type = 1; end
if ~isempty(nonzeros(M(:,maxy))) type = 2; end
if ~isempty(nonzeros(M(maxx,:))) type = 3; end
if ~isempty(nonzeros(M(:,1))) type = 4; end