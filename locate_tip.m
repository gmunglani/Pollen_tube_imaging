function [boundb, tip_final, tip_new, diam, maxy, center, phin, axes, stats, toln, major] = locate_tip(H,tol)

% Extract image boundary (longest boundary)
I = bwboundaries(H,'holes');
for x = 1:numel(I)
    tempbw(x) = size(I{x},1);
end
[pos posI] = max(tempbw);
bound = I{posI};

stats = regionprops(H,'Orientation','MajorAxisLength', 'BoundingBox', ...
    'MinorAxisLength', 'Eccentricity', 'Centroid','Area','FilledImage');

% Fit an ellipse to the entire image and get the maximum point
major = [stats.Centroid(2) + stats.MajorAxisLength*0.5 * sin(pi*stats.Orientation/180) ...
    stats.Centroid(1) - stats.MajorAxisLength*0.5 * cos(pi*stats.Orientation/180)];

% Remove points on the extreme right
maxy = size(H,2);
rem = find(bound(:,2) == maxy); bound(rem,:) = [];

% Find points on the convex hull
hullo = convhull(bound(:,1),bound(:,2));
ver = [bound(hullo,1) bound(hullo,2)];

% Find the diameter, midpoint at the cutoff along with the positions along
% the entire boundary
yedge = find(ver(:,2) == maxy-1);
[diam start edge] = edge_quant(ver,yedge);

% Tip finding algorithm
ybound = find(bound(:,2) == maxy-1);
boundc = circshift(bound,-ybound(1));
flag_tol = false; toln = tol;
while(flag_tol == 0)
    tip_new = [];
    for i = 1:length(bound)
        dist_val = pdist2(boundc(i,:),major(1,:));
        if (dist_val < toln*diam) tip_new = [tip_new; boundc(i,:)]; end
    end
    
    [tip_final,center,phin,axes,tip_check,fix, flag_tol] = ellipse_data(tip_new);
    cacl = axes(1) - axes(2);
    toln = toln - 0.1;
end

% Shift the entire boundary vector center at final tip
posxf = []; posxy = []; boundb = [];
posxf = find(boundc(:,1) == tip_final(1));
posyf = find(boundc(:,2) == tip_final(2));
interf = intersect(posxf,posyf);
boundb = circshift(boundc,(-ceil(length(boundc)*0.5)-interf(1)));
