% Gautam Munglani
% 15.08.2017

clear variables
close all

% Image path and input (MAKE SURE YOU HAVE THE RIGHT PATH)
fname_YFP ='J:\LRX_calcium\Data48_YFP.tif';
fname_CFP = 'J:\LRX_calcium\Data48_CFP.tif';

% Input parameters
thresh1 = 0.44; % Threshold value
pixelsize = 0.1; % Pixel to um conversion
gauss = 1.5; % Gauss filter option
npixel = 6; % Number of pixels difference from start point for straight line fit
stp = 1; % Start frame number
smp = 5; % End frame number

% Bleaching options
decay = -0.0008; % Bleaching decay (single exp decay, between -0.0005 and -0.0008)
time_step = 5; % Frame rate of movie

% Spline options
nint = 100; % Number of points to fit for spline
nbreaks = 5; % Number of spline regions

% ROI options
ROI_type = 1; % No ROI = 0; Moving ROI = 1; Stationary ROI = 2
split = 1; % Split ROI along center line
circle = 0; % Circle ROI as fraction of diameter
starti = 0; % Rectangle ROI Start length / no pixelsize means percentage as a fraction of length of tube
stopi = 3; % Rectangle/Circle ROI Stop length / no pixelsize means percentage as a fraction of length of tube

% Kymo and movie options
movie = '170622_YC_2_ratio7.avi'; % Make movie file if string is not empty
framerate = 10; % Video frame rate
Cmin = 0.5; % Min pixel value
Cmax = 1.5; % Max pixel value
nkymo = 5; % Number of pixels line width average for kymograph (even number) (0 means no kymo)

% Diameter options 
diamcutoff = 3; % distance from tip for first diameter calculations (um)

% Registration
register = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get info about image files
info1 = imfinfo(fname_YFP);
num_images1 = numel(info1);
info2 = imfinfo(fname_CFP);
num_images2 = numel(info2);

% Image Registration
[optimizer, metric] = imregconfig('multimodal');
for count = 1:smp-stp+1
    A1 = imread(fname_YFP, count, 'Info', info1); A1 = A1 ./ exp(decay*time_step*count);
    A2 = imread(fname_CFP, count, 'Info', info2);

    BO1(:,:,count) = imgaussfilt(mat2gray(A1),gauss); 
    BO2(:,:,count) = imgaussfilt(mat2gray(A2),gauss);  
end

% Crop to the max size
[tmp,posback] = imcrop(BO2(:,:,end));

for count = 1:smp-stp+1
    count
    B1 = BO1(:,:,count);
    BB1 = imcrop(B1,posback);
    B1(B1<thresh1) = 0; BB1(BB1>thresh1) = 0;
    B1 = B1 - sum(sum(BB1))/length(find(BB1>0));
    B1(B1<0) = 0; BF1(:,:,count) = B1; 
        
    B2 = BO2(:,:,count);    
    BB2 = imcrop(B2,posback);
    B2(B2<thresh1) = 0; BB2(BB2>thresh1) = 0;
    B2 = B2 - sum(sum(BB2))/length(find(BB2>0));
    B2(B2<0) = 0;
    
    if (register == 1) B2 = imregister(B2,B1,'similarity',optimizer,metric); end
    BF2(:,:,count) = B2;
end

[tmp,posfront] = imcrop(BF2(:,:,end));

if ~isempty(movie)
    [Cmin, Cmax, BN1, BN2, M] = video_processing(movie,stp,smp,BF1,BF2,posfront,framerate,time_step,Cmax,Cmin);
end

if (ROI_type > 0 || nkymo > 0 || diamcutoff > 0)
    for count = 1:smp-stp+1
        count
        if isempty(movie)
            % Orient image
            
            BN1 = imcrop(BF1(:,:,count),posfront);
            BN2 = imcrop(BF2(:,:,count),posfront);
            
            type = find_orient(BN1);
            if (type == 1) BN1 = imrotate(BN1,-90); BN2 = imrotate(BN2,-90);
            elseif (type == 3) BN1 = imrotate(BN1,90); BN2 = imrotate(BN2,90);
            elseif (type == 4) BN1 = imrotate(BN1,180); BN2 = imrotate(BN2,180);
            end
            
            % Ratio image and mask
            C = BN1./BN2;
            C(C==Inf) = 0;
            C(isnan(C)) = 0;
            
            M(:,:,count) = C;
        else
            C = M(:,:,count);
        end
 
        E = imbinarize(C,thresh1);
        
        % Fit ellipse continuously until major axis near tip is detected
        J = E;
        for a = 1:size(E,2)
            if(size(J,2) > min(find(sum(E,1)>1))+1.25.*sum(E(:,end))) J(:,end) = [];
            else break;
            end
            
            stats = regionprops(J,'Orientation','MajorAxisLength', 'BoundingBox', ...
                'MinorAxisLength', 'Eccentricity', 'Centroid','Area','FilledImage');
            
            % Fit an ellipse to the entire image and get the maximum point
            major(a,:) = [stats(1).Centroid(2) + stats(1).MajorAxisLength*0.5 * sin(pi*stats(1).Orientation/180) ...
                stats(1).Centroid(1) - stats(1).MajorAxisLength*0.5 * cos(pi*stats(1).Orientation/180)];
            
        end
        
        % Extract image boundary (longest boundary)
        I = bwboundaries(E,'holes');
        for x = 1:numel(I)
            temp1(x) = size(I{x},1);
        end
        [posI posI] = max(temp1);
        bound = I{posI};
        
        % Remove points on the extreme right
        maxy = size(E,2);
        rem = find(bound(:,2) == maxy); bound(rem,:) = [];
        
        % Find points on the convex hull
        hullo = convhull(bound(:,1),bound(:,2));
        ver = [bound(hullo,1) bound(hullo,2)];
        
        % Find the diameter, midpoint at the cutoff along with the positions along
        % the entire boundary
        yedge = find(ver(:,2) == maxy-1);
        [diam start edge] = edge_quant(ver,yedge);
        vers = circshift(ver,-start);
        
        % Locate which ellipse major axis point to use
        ell = abs(major(:,1) - mean(edge(:,1)));
        [peakdo,locel] = findpeaks(ell);
        if (isempty(locel)) majpt = length(ell);
        else majpt = locel(end);
        end
        
        % Find the points on the convex hull within a quarter of the diameter from the
        % ellipse major axis tip
        tip = []; tipsx = []; tipsy = []; tips = []; tipex = []; tipey = []; tipe = []; tipm = []; dist = [];
        for i = 1:length(vers)
            dist = pdist2(vers(i,:),major(majpt,:));
            if (dist < 0.6*diam) tip = [tip; vers(i,:)]; end
        end
        
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
        tip_final(count,:) = bounda(rangetip(1)+posid,:);
        
        % Shift the entire boundary vector center at final tip
        posxf = []; posxy = []; boundb = [];
        posxf = find(bounda(:,1) == tip_final(count,1));
        posyf = find(bounda(:,2) == tip_final(count,2));
        interf = intersect(posxf,posyf);
        boundb = circshift(bounda,(-ceil(length(bounda)*0.5)-interf(1)));
        
        % Find the curves along the sides of the tubes
        total1 = []; total2 = [];
        range1 = ceil(length(boundb)*0.5):length(boundb);
        dist1 = diag(pdist2(boundb(range1,:),(ones(length(range1),1))*tip_final(count,:)));
        postotal1 = find(dist1 > diam*0.75)+range1(1)-1;
        total1(:,:) = boundb(postotal1,:);
        
        range2 = ceil(length(boundb)*0.5)-1:-1:1;
        dist2 = diag(pdist2(boundb(range2,:),(ones(length(range2),1))*tip_final(count,:)));
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
        
        % Check for straight lines in Y near the tip
        linel1 = find((abs(total1(:,2) - total1(1,2)))<=npixel);
        if(length(linel1) > 5)
            posl1 = max(linel1);
            lenl1 = abs(total1(posl1,1) - total1(1,1));
        else
            posl1 = 1;
            lenl1 = 0;
        end
        
        linel2 = find((abs(total2(:,2) - total2(1,2)))<=6);
        if(length(linel2) > 5)
            posl2 = max(linel2);
            lenl2 = abs(total2(posl2,1) - total2(1,1));
        else
            posl2 = 1;
            lenl2 = 0;
        end
        
        % Divide edge array into straight lines and curves
        totall1 = total1(1:posl1,:);
        totall2 = total2(1:posl2,:);
        totals1 = total1(posl1:end,:);
        totals2 = total2(posl2:end,:);
        
        % Correct for angle at the right extreme of the image
        [totals1 totals2] = angle_correction(totals1,totals2,maxy,diam);
        
        % Find the lengths of individual segments of the side curves
        dlens1 = diag(pdist2(totals1(1:end-1,:),totals1(2:end,:)));
        dlens2 = diag(pdist2(totals2(1:end-1,:),totals2(2:end,:)));
        
        % find the number of points along the side curves
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
        fit1 = splinefit(totals1(:,2),totals1(:,1),breaks1,3);
        fit2 = splinefit(totals2(:,2),totals2(:,1),breaks2,3);
        
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
        distct = pdist2(tip_final(count,:),linec(1,:));
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
            diamf_avg(count) = sum(diamf)/length(diamf);
        end
        
        if (ROI_type > 0)
            % Find ROI from centerline distance using percentages or distance
            if (pixelsize == 0)
                percent = distc./distc(end);
                start_length = abs(percent*100 - starti); [startpos startpos] = min(start_length);
                stop_length = abs(percent*100 - stopi); [stoppos stoppos] = min(stop_length);
                lencp(count) = sum(lenc(startpos:stoppos));
                
                start_length = abs(percent*100*lencp(count)/lencp(stp) - starti); [startpos startpos] = min(start_length);
                stop_length = abs(percent*100*lencp(count)/lencp(stp) - stopi); [stoppos stoppos] = min(stop_length);
            else
                start_length = abs(distc*pixelsize + distct - starti); [startpos startpos] = min(start_length);
                stop_length = abs(distc*pixelsize + distct - stopi); [stoppos stoppos] = min(stop_length);
            end
            
            % Project ROI length onto the side curves
            [startc1,stopc1] = closest_bound(total1,xy1(:,2),xy1(:,1),startpos,stoppos);
            [startc2,stopc2] = closest_bound(total2,xy2(:,2),xy2(:,1),startpos,stoppos);
            
            % Create masks for rectangles and circles, and include whether they are
            % normal, split or stationary
            if (ROI_type ~= 2 | count == 1)
                if (circle == 0)
                    roi = vertcat(total1(startc1:stopc1,:), total2(stopc2:-1:startc2,:));
                    if (starti < distct) roi = vertcat(boundb(postotal2(1):postotal1(2),:),roi); end
                    F = poly2mask(roi(:,2),roi(:,1),size(E,1),size(E,2));
                else
                    mask = zeros(size(E,1),size(E,2));
                    if (stopi < distct)
                        npts = round(distct/(distc(1)*pixelsize));
                        xcir = linspace(tip_final(count,2),linec(1,2),npts);
                        ycir = linspace(tip_final(count,1),linec(1,1),npts);
                        distcir = diag(pdist2([tip_final(count,1)*ones(npts,1) tip_final(count,2)*ones(npts,1)], [ycir; xcir]'));
                        
                        if (pixelsize > 0)
                            stop_length = abs(distcir.*pixelsize - stopi);
                        else
                            stop_length = abs(distcir - stopi);
                        end
                        [stoppos stoppos] = min(stop_length);
                        roi = [ycir(stoppos) xcir(stoppos)];
                        mask(round(ycir(stoppos)),round(xcir(stoppos))) = 1;
                    else
                        roi = [round(linec(stoppos,1)) round(linec(stoppos,2))];
                        mask(round(linec(stoppos,1)),round(linec(stoppos,2))) = 1;
                    end
                    F = bwdist(mask) >= 0.5*circle.*diam;
                    F = imcomplement(F);
                end
                
                if (split == 1)
                    if (circle > 0)
                        stoppos = length(linec); stopc1 = length(total1); stopc2 = length(total2);
                    end
                    roi1 = vertcat(total1(startc1:stopc1,:), linec(stoppos:-1:startpos,:));
                    roi2 = vertcat(total2(startc2:stopc2,:), linec(stoppos:-1:startpos,:));
                    if (starti < distct)
                        roi1 = vertcat(boundb(range2(1):postotal1(2),:),roi1,boundb(range2(1),:));
                        roi2 = vertcat(boundb(range2(1):-1:postotal2(1),:),roi2,boundb(range2(1),:));
                    end
                    FS1 = F.*poly2mask(roi1(:,2),roi1(:,1),size(E,1),size(E,2));
                    FS2 = F.*poly2mask(roi2(:,2),roi2(:,1),size(E,1),size(E,2));
                end
            end
            % Create mask for intensity calculations
            maskim = cast(F, class(E));
            maskim_tmp = repmat(maskim, [1 1 3]);
            for w = 1:3
                G(:,:,w) = C .* maskim_tmp(:,:,w);
            end
            
            % Calculate average intensities and pixel numbers
            pixelnum(count) = sum(sum(F));
            intensity_avg(count) = sum(sum(C.*F))/pixelnum(count);
            intensityB1_avg(count) = sum(sum(BN1.*F))/pixelnum(count);
            intensityB2_avg(count) = sum(sum(BN2.*F))/pixelnum(count);
            
            % Repeat mask making and intensity averaging for the split case
            if (split == 1)
                maskim1 = cast(FS1, class(E));
                maskim1_tmp = repmat(maskim1, [1 1 3]);
                for w = 1:3
                    GS1(:,:,w) = C .*  maskim1_tmp(:,:,w);
                end
                
                pixelnumS1(count) = sum(sum(FS1));
                intensityS1_avg(count) = sum(sum(C.*FS1))/pixelnumS1(count);
                intensityB1S1_avg(count) = sum(sum(BN1.*FS1))/pixelnumS1(count);
                intensityB2S1_avg(count) = sum(sum(BN2.*FS1))/pixelnumS1(count);
                
                maskim2 = cast(FS2, class(E));
                maskim2_tmp = repmat(maskim2, [1 1 3]);
                for w = 1:3
                    GS2(:,:,w) = C .*  maskim2_tmp(:,:,w);
                end
                
                pixelnumS2(count) = sum(sum(FS2));
                intensityS2_avg(count) = sum(sum(C.*FS2))/pixelnumS2(count);
                intensityB1S2_avg(count) = sum(sum(BN1.*FS2))/pixelnumS2(count);
                intensityB2S2_avg(count) = sum(sum(BN2.*FS2))/pixelnumS2(count);
            end   
        end
        
        % Kymograph
        if (count == smp && nkymo > 0)
            % Find projection outside tip
            linect = vertcat(tip_final(smp,:),linec(1,:));
            fitct = polyfit(linect(:,2),linect(:,1),1);
            pointe(1,2) = 2*tip_final(smp,2)-linec(1,2);
            pointe(1,1) = fitct(1)*pointe(1,2) + fitct(2);
            linecte(:,:,1) = [vertcat(pointe(1,1), linec(:,1)), vertcat(pointe(1,2), linec(:,2))];
            
            % Average number of points across the tube width in kymo based on orientation
            for a = 2:nkymo
                ind = floor(a*0.5);
                if (start_nfitc(1) < 0)
                    if (mod(a,2) == 0) linecte(:,:,a) = [vertcat(pointe(1,1)+ind, linec(:,1)+ind), vertcat(pointe(1,2)-ind, linec(:,2)-ind)];
                    else linecte(:,:,a) = [vertcat(pointe(1,1)-ind, linec(:,1)-ind), vertcat(pointe(1,2)+ind, linec(:,2)+ind)];
                    end
                else
                    if (mod(a,2) == 0) linecte(:,:,a) = [vertcat(pointe(1,1)+ind, linec(:,1)+ind), vertcat(pointe(1,2)+ind, linec(:,2)+ind)];
                    else linecte(:,:,a) = [vertcat(pointe(1,1)-ind, linec(:,1)-ind), vertcat(pointe(1,2)-ind, linec(:,2)-ind)];
                    end
                end
            end
            for n = 1:count
                for a = 1:nkymo
                    L(:,:) = M(:,:,n)./Cmax;
                    L = (L-Cmin)./(1-Cmin);
                    L(L<0) = 0;
                    kymo(:,a) = improfile(L, linecte(:,2,a), linecte(:,1,a), round(sum(diag(pdist2(linecte(1:end-1,:,1),linecte(2:end,:,1))))));
                end
                kymo_avg(:,n) = mean(kymo(:,:),2);
            end
            kymo_avg(find(kymo_avg<0)) = 0;
        end
        
        % Images and mask visualization
        if (count == 1 || count == smp)
            figure
            subplot(2,3,1)
            imshow(B1);
            subplot(2,3,2)
            imshow(B2);
            subplot(2,3,3)
            imshow(C./max(M(:)));
            subplot(2,3,4)
            imshow(E);
            if (split == 1 && ROI_type > 0)
                subplot(2,3,5)
                imshow(FS1);
                subplot(2,3,6)
                imshow(FS2);
            elseif (ROI_type > 0)
                subplot(2,3,5)
                imshow(F);
            end
            
            % Plotting of data
            figure
            hold on
            plot(tip(:,2), tip(:,1), 'g*');
            plot(boundb(:,2), boundb(:,1), 'g', 'LineWidth', 2);
            plot(tip_final(count,2), tip_final(count,1), 'c+', 'LineWidth', 5);
            plot(major(majpt,2),major(majpt,1),'y+', 'LineWidth', 5);
            ellipse_view(stats);
            plot(total1(:,2), total1(:,1), 'k', 'LineWidth', 2);
            plot(total2(:,2), total2(:,1), 'b', 'LineWidth', 2);
            
            quiver(xc,yc,-dy,dx, 0, 'm')
            plot(xc, yc, 'c*', 'LineWidth', 0.5);
            
            circledraw([round(major(majpt,2)) round(major(majpt,1))],round(diam),100,':');
            if (ROI_type > 0)
                if (circle == 0)
                    plot(roi(:,2),roi(:,1),'r','LineWidth', 2);
                    if (split == 1)
                        plot(roi1(:,2),roi1(:,1),'r','LineWidth', 2);
                        plot(roi2(:,2),roi2(:,1),'r','LineWidth', 2);
                    end
                else
                    circledraw([roi(:,2) roi(:,1)],round(0.5*circle*diam),100,'r');
                end
            end
            
            if (diamcutoff > 0)                
                for i = 1:length(cut)
                    plot([xy1f(cut(i),2) xy2f(cut(i),2)],[xy1f(cut(i),1) xy2f(cut(i),1)],'k')
                    hold on
                end
                
                axis equal
                figure
                plot(distctf,diamf,'b')
            end
        end
    end
end
% Final tip movement/pixel number on a frame basis
if (diamcutoff > 0)
    figure
    subplot(1,2,1)
    plot(tip_final(:,2),tip_final(:,1))
    axis([0 size(E,2) 0 size(E,1)])
    title('Tip Final Position');

    subplot(1,2,2)
    plot(1:1:count,diamf_avg,'b')
    title('Average Diameter')
end

if (ROI_type > 0)
    figure    
    plot(1:1:count,pixelnum,'b')
    title('Pixel Number')
end

% Total intensity plots
if (ROI_type > 0)
    figure
    title('Average Intensities')
    if (split == 1) subplot(1,3,1); end
    hold on
    plot(1:1:count,intensity_avg,'g')
    plot(1:1:count,intensityB1_avg,'r')
    plot(1:1:count,intensityB2_avg,'b')
    
    if (split == 1)
        subplot(1,3,2)
        hold on
        plot(1:1:count,intensityS1_avg,'g')
        plot(1:1:count,intensityB1S1_avg,'r')
        plot(1:1:count,intensityB2S1_avg,'b')
        
        subplot(1,3,3)
        hold on
        plot(1:1:count,intensityS2_avg,'g')
        plot(1:1:count,intensityB1S2_avg,'r')
        plot(1:1:count,intensityB2S2_avg,'b')
    end
end

% Kymograph
if (nkymo > 0)
    figure
    map = colormap(jet(255));
    map = vertcat([0 0 0],map);
    imshow(uint8(kymo_avg.*255),map)
end