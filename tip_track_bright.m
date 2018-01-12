clear all
close all

% Image path and input
fname = '/home/gm/Documents/Work/Images/BF_tubes/rbohH33_adj.tiff';
info = imfinfo(fname);
num_images = numel(info);

% Input parameters
tol = 4; % Tolerance for tip finding algorithm (multiplier for circle diameter)
pixelsize = 0.1; % Pixel to um conversion
gauss = 2; % Gaussian smoothing

% Measurements
range_thresh = [-2,2]; % Multiplier for threshold
interval_thresh = 0.05; % Background threshold interval. 
npixel = 6; % Number of pixels difference from start point for straight line fit
hole_close = 1; % Hole closing algorithm
diamcutoff = 12; % Distance from tip for first diameter calculations (um)

% Spline options
nint = 100; % Number of points to fit for spline
nbreaks = 5; % Number of spline regions

% Reading the image, cropping and adding a blurring filter and thresholding
A = imread(fname);
B = imcrop(A);
C = imgaussfilt(mat2gray(B),gauss);
cutoff = multithresh(C,2);

e1 = strel('rectangle',[8 8]);
e2 = strel('disk',8);

h = figure;
img_count = 0;
for a = range_thresh(1):range_thresh(2)
    for b = range_thresh(1):range_thresh(2)
        range = [cutoff(1)+interval_thresh*a cutoff(2)-interval_thresh*b];
        D = imquantize(C,range);
        E = D; E(E==3) = 1; E(E==2) = 0;
        
        % Structuring elements and removing un-connected noise
        F = imdilate(E,e1); F = imerode(F,e2);
        G = bwareafilt(logical(F),1);
        
        % Orient image
        if (img_count == 0)
            type = find_orient(G);
            if (type == 1) G = imrotate(G,-90); C = imrotate(C,-90);
            elseif (type == 3) G = imrotate(G,90); C = imrotate(C,90);
            elseif (type == 4) G = imrotate(G,180); C = imrotate(C,180);
            end
        end
        
        Gmax = max(find(G(:,end)==1));
        Gmin = min(find(G(:,end)==1));
        H = imfill(drawline(G,Gmin,size(G,2),Gmax,size(G,2),1),'holes');
        
        img_count = img_count+1;
        thresh = range_thresh(1):1:range_thresh(2);
        len_thresh = length(thresh);
        subplot(len_thresh,len_thresh,img_count);
        imshowpair(C,H);
        if (mod(img_count-1,len_thresh) == 0) 
            posy = num2str(thresh(floor(img_count/len_thresh)+1));
            ylabel(posy); 
        end
        if ((img_count+len_thresh-1) >= len_thresh^2) 
            posx = num2str(thresh(img_count+len_thresh-len_thresh^2));
            xlabel(posx); 
        end
    end
end
mult_thresh1 = input(['Choose image with the best mask quality (row) [' num2str(range_thresh(1)) ' to ' num2str(range_thresh(2)) ']: ']);
mult_thresh2 = input(['Choose image with the best mask quality (column) [' num2str(range_thresh(1)) ' to ' num2str(range_thresh(2)) ']: ']);
close all

range = [cutoff(1)+interval_thresh*mult_thresh1 cutoff(2)-interval_thresh*mult_thresh2];
D = imquantize(C,range);
E = D; E(E==3) = 1; E(E==2) = 0;

% Structuring elements and removing un-connected noise
F = imdilate(E,e1); F = imerode(F,e2);
G = bwareafilt(logical(F),1);

% Orient image
type = find_orient(G);
if (type == 1) G = imrotate(G,-90); C = imrotate(C,-90);
elseif (type == 3) G = imrotate(G,90); C = imrotate(C,90);
elseif (type == 4) G = imrotate(G,180); C = imrotate(C,180);
end

Gmax = max(find(G(:,end)==1));
Gmin = min(find(G(:,end)==1));
H = imfill(drawline(G,Gmin,size(G,2),Gmax,size(G,2),1),'holes');
figure; imshowpair(C,H);

% Locate tip
[boundb, tip_final, tip_new, diam, maxy, center, phin, axes, stats, toln, major] = locate_tip(H,tol);

% Hole closing
if (hole_close == 1)
    val = pdist2(tip_final,tip_new);
    
    change = movmean(diff(val),5);
    for k = 1:length(change)
        if(change(k) > 0) turnin = k; break; end
    end
    for k = length(change):-1:1
        if(change(k) < 0) turnout = k+1; break; end
    end
    
    J = drawline(H,tip_new(turnin,1),tip_new(turnin,2),tip_new(turnout,1),tip_new(turnout,2),1);
    K = imfill(J,'holes');
    
    boundb = []; tip_new = [];
    [boundb, tip_final, tip_new, diam, maxy, center, phin, axes, stats, toln, major] = locate_tip(K,tol);
end 

% Find the curves along the sides of the tubes
total1 = []; total2 = [];
range1 = ceil(length(boundb)*0.5):length(boundb);
dist1 = diag(pdist2(boundb(range1,:),(ones(length(range1),1))*tip_final));
postotal1 = find(dist1 > diam*0.75)+range1(1)-1;
total1(:,:) = boundb(linspace(postotal1(1),postotal1(end),length(postotal1)),:);

range2 = ceil(length(boundb)*0.5)-1:-1:1;
dist2 = diag(pdist2(boundb(range2,:),(ones(length(range2),1))*tip_final));
postotal2 = range2(1)-find(dist2 > diam*0.75)+1;
total2(:,:) = boundb(linspace(postotal2(1),postotal2(end),length(postotal2)),:);

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
%plot(major(:,2),major(:,1),'k*', 'LineWidth', 4);
ellipse_view(center,phin,axes);
%ellipse_view(stats.Centroid(2:-1:1),pi*stats.Orientation/180,[stats.MajorAxisLength/2 stats.MinorAxisLength/2]);
circledraw([round(major(:,2)) round(major(:,1))],round(toln*diam),100,':');
plot(tip_final(2), tip_final(1), 'r*', 'LineWidth', 4);
plot(total1(:,2), total1(:,1), 'k', 'LineWidth', 2);
plot(total2(:,2), total2(:,1), 'b', 'LineWidth', 2);
plot(tip_new(:,2), tip_new(:,1), 'm', 'LineWidth', 2);
quiver(xc,yc,-dy,dx, 0, 'm')
plot(xc, yc, 'c*', 'LineWidth', 0.5);
%plot(center(:,2),center(:,1),'k*','LineWidth',2)
set(gca,'Ydir','reverse')

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

Average_diam = mean(diamf)
