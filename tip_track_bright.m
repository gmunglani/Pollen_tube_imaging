clear all
close all

% Image path and input
fname = '';
info = imfinfo(fname);
num_images = numel(info);

% Input parameters
tol = 0.6; % Tolerance for tip finding algorithm (multiplier for circle diameter)
pixelsize = 0.1; % Pixel to um conversion
gauss = 2;
mult_thresh = 1.14;
npixel = 6; % Number of pixels difference from start point for straight line fit
diamcutoff = 4; % Distance from tip for first diameter calculations (um)

% Spline options
nint = 100; % Number of points to fit for spline
nbreaks = 5; % Number of spline regions

% Reading the image, cropping and adding a blurring filter and thresholding
A = imread(fname);
B = imcrop(A);
C = imgaussfilt(mat2gray(B),gauss);
cutoff = multithresh(C,2);
D = imquantize(C,[cutoff(1)*mult_thresh cutoff(2)/mult_thresh]);
E = D; E(E==3) = 1; E(E==2) = 0;

% Structuring elements and removing un-connected noise
e1 = strel('rectangle',[8 8]);
e2 = strel('disk',8);
F = imdilate(E,e1); F = imerode(F,e2);
G = imfill(F,'holes');
H = bwareafilt(logical(G),1); 

% Orient image
type = find_orient(H);
if (type == 1) H = imrotate(H,-90); 
elseif (type == 3) H = imrotate(H,90); 
elseif (type == 4) H = imrotate(H,180); 
end
        
% Extract image boundary (longest boundary)
I = bwboundaries(H,'holes');
for x = 1:numel(I)
    tempbw(x) = size(I{x},1);
end
[pos posI] = max(tempbw);
bound = I{posI};

stats_orig = regionprops(H,'Orientation','MajorAxisLength', 'BoundingBox', ...
    'MinorAxisLength', 'Eccentricity', 'Centroid','Area','FilledImage');

Mhmax = max(find(H(:,end)==1));
Mhmin = min(find(H(:,end)==1));
Mwmax = min(find(H(floor(stats_orig.BoundingBox(2)+stats_orig.BoundingBox(4)),:)==1));
Mwmin = min(find(H(ceil(stats_orig.BoundingBox(2)),:)==1));

Mangle_up = atan2(abs(stats_orig.BoundingBox(2)+stats_orig.BoundingBox(4)-Mhmax),size(H,2)-Mwmax)*180/pi;
Mangle_do = atan2(abs(Mhmin-stats_orig.BoundingBox(2)),size(H,2)-Mwmin)*180/pi;

Mangle = max(abs(Mangle_up),abs(Mangle_do));
Mratio = stats_orig.BoundingBox(4)/(Mhmax-Mhmin);
Mxdist = 1.55 - min(0.3*(Mratio)*Mangle/15,0.25);

J = H(:,1:floor(stats_orig.BoundingBox(1)+Mxdist*sum(H(:,end))));

% Filled in image with the given outline
stats = regionprops(J,'Orientation','MajorAxisLength', 'BoundingBox', ...
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
vers = circshift(ver,-start);

% Find the points on the convex hull within a ratio of the diameter from the
% ellipse major axis tip
tip = []; tipsx = []; tipsy = []; tips = []; tipex = []; tipey = []; tipe = []; tipm = []; dist = [];
for i = 1:length(vers)
    dist = pdist2(vers(i,:),major(1,:));
    if (dist < tol*diam) tip = [tip; vers(i,:)]; end
end
if isempty(tip) error('Increase tolerance'); end

% Shift array to start at first point of tip area
% Get the start(s), end(e) and mid(m) of the tip area
tipsx = [find(bound(:,1) == tip(1,1))];
tipsy = [find(bound(:,2) == tip(1,2))];
intertips = intersect(tipsx,tipsy);
bounda = circshift(bound,-intertips(1));
tips = bounda(1,:);

tipex = [find(bounda(:,1) == tip(end,1))];
tipey = [find(bounda(:,2) == tip(end,2))];
intertipe = intersect(tipex,tipey);
tipe = bounda(intertipe(1),:);

tipm = bounda(round((length(bounda)+intertipe)/2),:);
rangetip = intertipe(1):length(bounda);

% Calculate the differential of all points in the tip region from the
% start, mid and end points
[dlenm] = dist_smooth(tipm(1,:),rangetip,bounda,1,5);
[dlene] = dist_smooth(tipe(1,:),rangetip,bounda,1,5);
[dlens] = dist_smooth(tips(1,:),rangetip,bounda,1,5);

% Locate the final tip by analysing the data
posid = round(length(dlenm)*0.5);
for k = 1:length(dlenm)
    if((dlene(k) + dlens(k)) < 0) posid = k; break; end
end
tip_final = bounda(rangetip(1)+posid,:);

% Shift the entire boundary vector center at final tip
posxf = []; posxy = []; boundb = [];
posxf = find(bounda(:,1) == tip_final(1));
posyf = find(bounda(:,2) == tip_final(2));
interf = intersect(posxf,posyf);
boundb = circshift(bounda,(-ceil(length(bounda)*0.5)-interf(1)));

% Find the curves along the sides of the tubes
total1 = []; total2 = [];
range1 = ceil(length(boundb)*0.5):length(boundb);
dist1 = diag(pdist2(boundb(range1,:),(ones(length(range1),1))*tip_final));
postotal1 = find(dist1 > diam*0.75)+range1(1)-1;
total1(:,:) = boundb(postotal1,:);

range2 = ceil(length(boundb)*0.5)-1:-1:1;
dist2 = diag(pdist2(boundb(range2,:),(ones(length(range2),1))*tip_final));
postotal2 = range2(1)-find(dist2 > diam*0.75)+1;
total2(:,:) = boundb(postotal2,:);

% Ensure that both curves reach maxy
if (max(total1(:,2)) < (maxy-1))
    while(max(total1(:,2)) < (maxy-1))
        total1 = vertcat(total1,total2(end,:));
        total2(end,:) = [];
    end
elseif (max(total2(:,2)) < (maxy-1))
    while(max(total2(:,2)) < (maxy-1))
        total2 = vertcat(total2, total1(end,:));
        total1(end,:) = [];
    end
end

if(abs(total1(end,1) - total2(end,1)) < 0.75*diam)
    total1(find(total1(:,2) >= max(total1(:,2))),:) = [];
    total2(find(total2(:,2) >= max(total2(:,2))),:) = [];
end

% Check for straight lines in Y near thtip_finale tip
linel1 = find((abs(total1(:,2) - total1(1,2)))<=npixel);
if(length(linel1) > 5)
    posl1 = max(linel1);
    lenl1 = abs(total1(posl1,1) - total1(1,1));
else
    posl1 = 1;
    lenl1 = 0;
end

linel2 = find((abs(total2(:,2) - total2(1,2)))<=npixel);
if(length(linel2) > 5)
    posl2 = max(linel2);
    lenl2 = abs(total2(posl2,1) - total2(1,1));
else
    posl2 = 1;
    lenl2 = 0;
end

% Divide edge array into straight linestip_final and curves
totall1 = total1(1:posl1,:);
totall2 = total2(1:posl2,:);
totals1 = total1(posl1:end,:);
totals2 = total2(posl2:end,:);

% Correct for angle at the right extreme of the image
[totals1 totals2] = angle_correction(totals1,totals2,maxy,diam);

% Find the lengths of individual segments of the side curves
dlens1 = diag(pdist2(totals1(1:end-1,:),totals1(2:end,:)));
dlens2 = diag(pdist2(totals2(1:end-1,:),totals2(2:end,:)));

% find the number of points along the side curvestip_final
nints1 = floor(nint*sum(dlens1)/(sum(dlens1)+lenl1));
nints2 = floor(nint*sum(dlens2)/(sum(dlens2)+lenl2));

% Finding the breaks for the spline fit
breaks1(1) = totals1(1,2); breaks2(1) = totals2(1,2);
for i = 1:nbreaks
    for j = 1:length(totals1)-1
        if (sum(dlens1(1:j))./sum(dlens1) >= i/nbreaks) breaks1(i+1) = totals1(j+1,2); break;end
    end
    for j = 1:length(totals2)-1
        if (sum(dlens2(1:j))./sum(dlens2) >= i/nbreaks) breaks2(i+1) = totals2(j+1,2); break;end
    end
end

% Fit splines to side curves
fit1 = splinefit(totals1(:,2),totals1(:,1),breaks1);
fit2 = splinefit(totals2(:,2),totals2(:,1),breaks2);

% Find x positions for fits and evaluate the spline at the given locations
xx1 = []; xx2 = []; yy1 = []; yy2 = []; distc = []; distct = []; distcf = []; diamf = [];
x1 = linspace(0,sum(dlens1),nints1);
for i = 1:length(x1)
    for j = 1:length(totals1)-1
        if (sum(dlens1(1:j)) >= x1(i)) xx1(i) = totals1(j+1,2); break;end
    end
end
%  xx1 = space_out(xx1);
yy1 = ppval(fit1,xx1);
yy1l = linspace(totall1(1,1),totall1(end,1),nint-nints1+1);
xx1l = linspace(totall1(1,2),totall1(end,2),nint-nints1+1);
yy1 = [yy1l(1:end-1) yy1];
xx1 = [xx1l(1:end-1) xx1];

x2 = linspace(0,sum(dlens2),nints2);
for i = 1:length(x2)
    for j = 1:length(totals2)-1
        if (sum(dlens2(1:j)) >= x2(i)) xx2(i) = totals2(j+1,2); break;end
    end
end
% xx2 = space_out(xx2);
yy2 = ppval(fit2,xx2);
yy2l = linspace(totall2(1,1),totall2(end,1),nint-nints2+1);
xx2l = linspace(totall2(1,2),totall2(end,2),nint-nints2+1);
yy2 = [yy2l(1:end-1) yy2];
xx2 = [xx2l(1:end-1) xx2];

% Average splines to get the center line of the tube
xc = 0.5*(xx1 + xx2);
yc = 0.5*(yy1 + yy2);
linec = [yc' xc'];

% Find the length of the center line including and excluding the tip
distct = pdist2(tip_final,linec(1,:));
if (pixelsize > 0) distct = pixelsize*distct; end

lenc = diag(pdist2(linec(1:end-1,:),linec(2:end,:)));
for s = 1:length(lenc) distc(s) = sum(lenc(1:s));end

% Calculate the gradient of the center line to get the normals
dx = gradient(xc);
dy = gradient(yc);

% Finding the points where the normals hit the edge curves
cross1 = 1; cross2 = 1; poscross1 = []; poscross2 = [];
for n = 1:nint
    nfitc = polyfit(vertcat(xc(n),(xc(n) - dy(n))),vertcat(yc(n),(yc(n) + dx(n))),1);
    if (n == 1) start_nfitc(:,:) = nfitc(:,:); end
    
    edge1 = total1(:,1) - nfitc(1).*total1(:,2) - nfitc(2);
    [cross1 cross1] = min(abs(edge1));
    poscross1(n) = cross1;
    
    edge2 = total2(:,1) - nfitc(1).*total2(:,2) - nfitc(2);
    [cross2 cross2] = min(abs(edge2));
    poscross2(n) = cross2;
end


% Ensure that all overlapping diameter lines are shifted backwards to
% ensure continuity
xy1 = []; xy2 = []; xy1f = []; xy2f = [];
xy1(:,:) = total1(poscross1,:);
xy2(:,:) = total2(poscross2,:);

[poscross1, poscross2, distcf] = line_continuity(poscross1,poscross2,1,distc);
[poscross1, poscross2, distcf] = line_continuity(poscross1,poscross2,2,distcf);

xy1f(:,:) = total1(poscross1,:);
xy2f(:,:) = total2(poscross2,:);

if (diamcutoff >0)
    % Find actual diameter
    if (pixelsize > 0) diamoff = find((pixelsize*distcf) > diamcutoff-distct); fcut = diamoff(1);
    else fcut = floor(size(xy1f,1)*0.05);
    end
    
    cut = fcut:ceil(size(xy1f,1)*0.95);
    diamf = diag(pdist2(xy1f(cut,:),xy2f(cut,:)));
    if (pixelsize > 0) diamf = pixelsize*diamf; distcf = pixelsize*distcf; end
    distctf = distcf(cut) + distct;
    diamf_avg = sum(diamf)/length(diamf);
end

figure
hold on
plot(boundb(:,2), boundb(:,1), 'g', 'LineWidth', 2);
plot(major(:,2),major(:,1),'k*', 'LineWidth', 4);
ellipse_view(stats);
circledraw([round(major(:,2)) round(major(:,1))],round(tol*diam),100,':');
plot(tip_final(2), tip_final(1), 'r*', 'LineWidth', 4);
plot(tip(:,2), tip(:,1), 'm', 'LineWidth', 2);
plot(total1(:,2), total1(:,1), 'k', 'LineWidth', 2);
plot(total2(:,2), total2(:,1), 'b', 'LineWidth', 2);
quiver(xc,yc,-dy,dx, 0, 'm')
plot(xc, yc, 'c*', 'LineWidth', 0.5);

if (diamcutoff > 0)
    for i = 1:length(cut)
        plot([xy1f(cut(i),2) xy2f(cut(i),2)],[xy1f(cut(i),1) xy2f(cut(i),1)],'k')
        hold on
    end
    axis equal

    figure
    plot(distctf,diamf,'b')
    title('Diameter with distance from the tip');
end

