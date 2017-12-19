clear all
close all

% Image path and input (MAKE SURE YOU HAVE THE RIGHT PATH)
fname_YFP = '/home/gm/Documents/Work/Images/Ratio_tubes/YC11_YFP.tif';
fname_CFP = '/home/gm/Documents/Work/Images/Ratio_tubes/YC11_CFP.tif';

% Input parameters
tol = 1; % Tolerance for tip finding algorithm (multiplier for circle diameter)
analysis = 1; % Turn on analysis mode
details = 0;  % Show histograms of results in the end
pixelsize = 0.1; % Pixel to um conversion
npixel = 6; % Number of pixels difference from start point for straight line fit
stp = 1; % Start frame number
smp = 5; % End frame number

% Spline options
nint = 100; % Number of points to fit for spline
nbreaks = 5;% Number of spline regions

% ROI options
ROItype = 1; % No ROI = 0; Moving ROI = 1; Stationary ROI = 2
split = 0; % Split ROI along center line
circle = 0; % Circle ROI as fraction of diameter
starti = 0; % Rectangle ROI Start length / no pixelsize means percentage as a fraction of length of tube
stopi = 30; % Rectangle/Circle ROI Stop length / no pixelsize means percentage as a fraction of length of tube

% Kymo and movie options
path = '~/Documents/Scripts/MATLAB/tip_results/test8'; % Make movie file if string is not empty
video_intensity = 0; % Video intensity
video_plot = 1; % Video of tip detection
timestep = 0.25; % Frame rate of movie
Cmin = 2.1; % Min pixel value
Cmax = 2.7; % Max pixel value
nkymo = 0; % Number of pixels line width average for kymograph (even number) (0 means no kymo)

% Diameter options 
diamcutoff = 0; % Distance from tip for first diameter calculations (um)

% Registration
register = 1; % Register image
union = 1; % Take the union of the two image masks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frame range
if  exist([path '.mat'],'file') == 2    
    load([path '.mat']);
    analysis = 0;
    disp('Analysis file exists');
else
    % Get info about image files
    info1 = imfinfo(fname_YFP);
    num_images1 = numel(info1);
    info2 = imfinfo(fname_CFP);
    num_images2 = numel(info2);
    
    for proc = stp:smp
        
    end
    
    % Crop region on the last frame
    Al = imread(fname_CFP, smp, 'Info', info1);
    Bl = imgaussfilt(mat2gray(Al),gauss);
    [tmp,posfront] = imcrop(Bl);
    [optimizer, metric] = imregconfig('multimodal');

    Bedge = zeros(1,smp); Bsum = zeros(1,smp);
    
    for count = stp:smp
        % Read image and add bleach correction
        disp(['Pre Processing:' num2str(count)]);
        
        A1 = imread(fname_YFP,count); 
        A2 = imread(fname_CFP,count);
        
        AF1 = imcrop(A1,posfront)
        AF2 = imcrop(A2,posfront)
        
        if (register == 1) BF2 = imregister(BF2,BF1,'translation',optimizer,metric); end
        
        % Orient image
        if (count==stp) type = find_orient(BF1); end
        if (type == 1) BN1 = imrotate(BN1,-90); BN2 = imrotate(BN2,-90);
        elseif (type == 3) BN1 = imrotate(BN1,90); BN2 = imrotate(BN2,90);
        elseif (type == 4) BN1 = imrotate(BN1,180); BN2 = imrotate(BN2,180);
        end
        
        if (count == stp)
            BT1=zeros(size(BN1,1),size(BN1,2),smp);
            BT2=zeros(size(BN2,1),size(BN2,2),smp);
        end
        
        % Union of images to ensure perfect overlap
        B = BN1.*BN2; B(B>0) = 1;
        while(sum(B(:,end-Bedge(count))) == 0)
            Bedge(count) = Bedge(count) + 1;
        end

        if (union == 1)
            BT1(:,:,count) = BN1.*B;
            BT2(:,:,count) = BN2.*B;
        else
            BT1(:,:,count) = BN1;
            BT2(:,:,count) = BN2;
        end
        
        Bsum(count) = sum(B(:));
        BN1sum(count) = nnz(BN1);
        BN2sum(count) = nnz(BN2);
        BF1sum(count) = nnz(BF1);
        BF2sum(count) = nnz(BF2);
    end
    
    BT1 = BT1(:,1:end-max(Bedge),:);
    BT2 = BT2(:,1:end-max(Bedge),:);
    
    BT1max = max(BT1(:));
    BT1 = BT1./BT1max;
    BT2 = BT2./BT1max; 
        
    % Create ratio image and output values
    M = BT1./BT2; % Set M is equal to the channel of interest
    M(M==Inf) = 0;
    M(isnan(M)) = 0;
end

if (analysis == 1)
    h = figure;
    figure(h);
    subplot(2,2,1)
    plot(stp:smp,Bsum,'b');
    hold on
    plot(stp:smp,BN1sum,'r');
    plot(stp:smp,BN2sum,'g');
    axis([0 smp 0.5*max(Bsum) max(Bsum)])
    title('Total Number of Pixels');

    [Mmin Mmax Mmin_prc Mmax_prc] = channel_analysis(M,smp);
    [B1min B1max B1min_prc B1max_prc] = channel_analysis(BT1,smp);
    [B2min B2max B2min_prc B2max_prc] = channel_analysis(BT2,smp);
    
    subplot(2,2,2)
    hold on
    plot(stp:smp,Mmax,'b*')
    plot(stp:smp,Mmax_prc,'bd')
    plot(stp:smp,Mmin,'r*')
    plot(stp:smp,Mmin_prc,'rd')
    grid on
    axis([0 smp 0 max(Mmax)])
    set(gca,'YMinorTick','on')
    title('Percentiles of Ratio image');

    subplot(2,2,3)
    hold on
    plot(stp:smp,B1max,'b*')
    plot(stp:smp,B1max_prc,'bd')
    plot(stp:smp,B1min,'r*')
    plot(stp:smp,B1min_prc,'rd')
    grid on      
    axis([0 smp 0 1])
    set(gca,'YMinorTick','on');
    title('Percentiles of YFP image');
    
    subplot(2,2,4)
    hold on
    plot(stp:smp,B2max,'b*')
    plot(stp:smp,B2max_prc,'bd')
    plot(stp:smp,B2min,'r*')
    plot(stp:smp,B2min_prc,'rd')
    grid on
    axis([0 smp 0 1])
    set(gca,'YMinorTick','on');
    title('Percentiles of CFP image');
    hold off
    save([path '.mat'],'BT1','BT2','M');
    savefig(h,path);
    error('Check pixel counts and min/max intensity');
end

% Make a movie and output min and max intensities of the whole stack
if (video_intensity == 1)
   video_processing(movie,stp,smp,BT1,BT2,framerate,timestep,Cmax,Cmin,M);
end

% Video plot
if (video_plot == 1)
    path = [path '_ellipse.avi'];
    V = VideoWriter(path);
    V.FrameRate = 1/timestep;
    open(V);
end

if (ROItype > 0 || nkymo > 0 || diamcutoff > 0)
    for count = stp:smp
        disp(['Image Analysis:' num2str(count)]);
        % Ratio image binarization
        C = M(:,:,count);
        D = bwareafilt(imbinarize(C,0.5),1);
        
        Ehmax = max(find(D(:,end)==1));
        Ehmin = min(find(D(:,end)==1));
        
        E = imfill(drawline(D,Ehmin,size(D,2),Ehmax,size(D,2),1),'holes');

        % Extract image boundary (longest boundary)
        I = bwboundaries(E,'holes');
        temp1 = [];
        for x = 1:numel(I)
            temp1(x) = size(I{x},1);
        end
        [posI posI] = max(temp1);
        bound = I{posI};
                
        stats_orig = regionprops(E,'Orientation','MajorAxisLength', 'BoundingBox', ...
                'MinorAxisLength', 'Eccentricity', 'Centroid','Area','FilledImage');
            
        Ewmax = min(find(E(floor(stats_orig.BoundingBox(2)+stats_orig.BoundingBox(4)),:)==1));
        Ewmin = min(find(E(ceil(stats_orig.BoundingBox(2)),:)==1));
        
        Eangle_up = atan2(abs(stats_orig.BoundingBox(2)+stats_orig.BoundingBox(4)-Ehmax),size(E,2)-Ewmax)*180/pi;    
        Eangle_do = atan2(abs(Ehmin-stats_orig.BoundingBox(2)),size(E,2)-Ewmin)*180/pi;    
        
        Eangle = max(abs(Eangle_up),abs(Eangle_do));
        Eratio = stats_orig.BoundingBox(4)/(Ehmax-Ehmin);
        Exdist = 1.5 - min(0.3*(Eratio)*Eangle/15,0.3); 
        
        J = E(:,1:floor(stats_orig.BoundingBox(1)+Exdist*sum(E(:,end))));
        stats = regionprops(J,'Orientation','MajorAxisLength', 'BoundingBox', ...
                'MinorAxisLength', 'Eccentricity', 'Centroid','Area','FilledImage');
        
        major(1,:) = [stats(1).Centroid(2) + stats(1).MajorAxisLength*0.5 * sin(pi*stats(1).Orientation/180) ...
                      stats(1).Centroid(1) - stats(1).MinorAxisLength*0.5 * cos(pi*stats(1).Orientation/180)];             
        
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
        if (count == stp) diamo = diam; end

        % Tip finding algorithm
        ybound = find(bound(:,2) == maxy-1);
        boundc = circshift(bound,-ybound(1));
        flag_tol = false; toln = tol;
        while(flag_tol == 0)
            tip_new = [];
            for i = 1:length(bound)
                dist_val = pdist2(boundc(i,:),major(1,:));
                if (dist_val < toln*diamo) tip_new = [tip_new; boundc(i,:)]; end
            end
       
            [tip_final(count,:),center,phin(count),axes,tip_check,fix(count), flag_tol] = ellipse_data(tip_new);
            cacl(count) = axes(1) - axes(2);
            toln = toln - 0.1;
        end

        % Shift the entire boundary vector center at final tip
        posxf = []; posxy = []; boundb = [];
        posxf = find(boundc(:,1) == tip_final(count,1));
        posyf = find(boundc(:,2) == tip_final(count,2));
        interf = intersect(posxf,posyf);
        boundb = circshift(boundc,(-ceil(length(boundc)*0.5)-interf(1)));
        
        % Find the curves along the sides of the tubes
        total1 = []; total2 = [];
        range1 = ceil(length(boundb)*0.5):length(boundb);
        dist1 = diag(pdist2(boundb(range1,:),(ones(length(range1),1))*tip_final(count,:)));
        postotal1 = find(dist1 > diamo*0.75)+range1(1)-1;
        total1(:,:) = boundb(postotal1,:);
        
        range2 = ceil(length(boundb)*0.5)-1:-1:1;
        dist2 = diag(pdist2(boundb(range2,:),(ones(length(range2),1))*tip_final(count,:)));
        postotal2 = range2(1)-find(dist2 > diamo*0.75)+1;
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
        
        % Check for straight lines in Y near the tip
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
        
        % Divide edge array into straight lines and curves
        totall1 = total1(1:posl1,:);
        totall2 = total2(1:posl2,:);
        totals1 = total1(posl1:end,:);
        totals2 = total2(posl2:end,:);
        
        % Correct for angle at the right extreme of the image
        [totals1 totals2] = angle_correction(totals1,totals2,maxy,diamo);
        
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
        
        % Remove first few points
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
        
        if (ROItype > 0)
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
            if (ROItype ~= 2 | count == stp)
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
                    F = bwdist(mask) >= 0.5*circle.*diamo;
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
            
            % Calculate average intensities and pixel numbers
            pixelnum(count) = sum(sum(F));
            intensity_tavg(count) = sum(sum(C))/nnz(C);
            intensity_avg(count) = sum(sum(C.*F))/pixelnum(count);
            intensityB1_avg(count) = sum(sum(BT1(:,:,count).*F))/pixelnum(count);
            intensityB2_avg(count) = sum(sum(BT2(:,:,count).*F))/pixelnum(count);
        end
        
        % Kymograph
        if (count == smp && nkymo > 0)
            % Find projection outside tip
            linect = vertcat(tip_final(count,:),linec(1,:));
            fitct = polyfit(linect(:,2),linect(:,1),1);
            pointe(1,2) = 2*tip_final(count,2)-linec(1,2);
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
                    L = (L-(Cmin/Cmax))./(1-(Cmin/Cmax));
                    L(L<0) = 0;
                    kymo(:,a) = improfile(L, linecte(:,2,a), linecte(:,1,a), round(sum(diag(pdist2(linecte(1:end-1,:,1),linecte(2:end,:,1))))));
                end
                kymo_avg(:,n) = mean(kymo(:,:),2);
            end
            kymo_avg(find(kymo_avg<0)) = 0;
        end
        
        if (video_plot == 1)
            % Plotting video of data
            h = figure;
            figure(h);
            hold on
            plot(boundb(:,2), boundb(:,1), 'g', 'LineWidth', 2);
            ellipse_view(center,phin(count),axes);
            circledraw([round(major(:,2)) round(major(:,1))],round((toln+0.1)*diamo),100,'k:');
            plot(tip_final(count,2), tip_final(count,1), 'm*', 'LineWidth', 4);
            plot(tip_check(:,2), tip_check(:,1), 'b*', 'LineWidth', 4);            
            
            plot(center(2),center(1),'b*', 'LineWidth',4);
            plot(total1(:,2), total1(:,1), 'k', 'LineWidth', 2);
            plot(total2(:,2), total2(:,1), 'b', 'LineWidth', 2);
            plot(tip_new(:,2), tip_new(:,1), 'm', 'LineWidth', 2);
            quiver(xc,yc,-dy,dx, 0, 'm')
            plot(xc, yc, 'c*', 'LineWidth', 0.5);

%             if (ROItype > 0)
%                 if (circle == 0)
%                     plot(roi(:,2),roi(:,1),'r','LineWidth', 2);
%                     if (split == 1)
%                         plot(roi1(:,2),roi1(:,1),'r','LineWidth', 2);
%                         plot(roi2(:,2),roi2(:,1),'r','LineWidth', 2);
%                     end
%                 else
%                     circledraw([roi(:,2) roi(:,1)],round(0.5*circle*diamo),100,'r');
%                 end
%             end
        
            axis([0 max(size(E,2),size(E,1)) 0 max(size(E,2),size(E,1))])
            title(['Frame:' num2str(count)])
            frame = getframe(gcf);
            writeVideo(V,frame);
            close(h);
        end
        
        if (diamcutoff > 0)
            for i = 1:length(cut)
                plot([xy1f(cut(i),2) xy2f(cut(i),2)],[xy1f(cut(i),1) xy2f(cut(i),1)],'k')
            end
            axis equal
            
            figure
            plot(distctf,diamf,'b')
            title('Diameter with distance from the tip');
        end
        

        
        % Images and mask visualization
        if (count == stp || count == smp)
            if(count == stp && details == 1) 
                Chist1 = reshape(C,1,size(C,1)*size(C,2)); 
                B1hist1 = reshape(BT1(:,:,count),1,size(BT1(:,:,count),1)*size(BT1(:,:,count),2));
                B2hist1 = reshape(BT2(:,:,count),1,size(BT2(:,:,count),1)*size(BT2(:,:,count),2));
                
                ChistF1 = reshape(C.*F,1,size(C,1)*size(C,2)); 
                B1histF1 = reshape(BT1(:,:,count).*F,1,size(BT1(:,:,count),1)*size(BT1(:,:,count),2));
                B2histF1 = reshape(BT2(:,:,count).*F,1,size(BT2(:,:,count),1)*size(BT2(:,:,count),2));
            else 
                Chist2 = reshape(C,1,size(C,1)*size(C,2)); 
                B1hist2 = reshape(BT1(:,:,count),1,size(BT1(:,:,count),1)*size(BT1(:,:,count),2));
                B2hist2 = reshape(BT2(:,:,count),1,size(BT2(:,:,count),1)*size(BT2(:,:,count),2));
                
                ChistF2 = reshape(C.*F,1,size(C,1)*size(C,2)); 
                B1histF2 = reshape(BT1(:,:,count).*F,1,size(BT1(:,:,count),1)*size(BT1(:,:,count),2));
                B2histF2 = reshape(BT2(:,:,count).*F,1,size(BT2(:,:,count),1)*size(BT2(:,:,count),2));
            end
            
            figure
            subplot(2,3,1)
            imshow(BT1(:,:,count));
            subplot(2,3,2)
            imshow(BT2(:,:,count));
            subplot(2,3,3)
            imshow(C./max(M(:)));
            subplot(2,3,4)
            imshow(E);
            if (split == 1 && ROItype > 0)
                subplot(2,3,5)
                imshow(FS1);
                subplot(2,3,6)
                imshow(FS2);
            elseif (ROItype > 0)
                subplot(2,3,5)
                imshow(F);
            end            
        end
    end
end

if (video_plot == 1) close(V); end

%Final tip movement/pixel number on a frame basis
if (diamcutoff > 0)
    figure
    subplot(1,2,1)
    plot(tip_final(:,2),tip_final(:,1))
    axis([0 size(E,2) 0 size(E,1)])
    set(gca, 'FontSize', 16)
    title('Tip Final Position');

    subplot(1,2,2)
    plot(1:length(diamf_avg),diamf_avg,'b')
    xlabel('Frame number', 'FontSize',12)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16)
    ylabel('Average Diameter', 'FontSize',12)
    yt = get(gca, 'YTick');
    title('Average Diameter')
    axis([0 max(count) 0 max(diamf_avg)])
end

if (ROItype > 0)
    figure    
    plot(1:length(pixelnum),pixelnum,'b')
    title('Pixel Number')
    axis([0 max(count) 0 max(pixelnum)+max(pixelnum)/10])
end

% Total intensity plots
if (ROItype > 0)
    figure
    if (split == 1) 
        subplot(1,3,2)
        hold on
        plot(1:length(intensityS1_avg),intensityS1_avg,'g')
        plot(1:length(intensityB1S1_avg),intensityB1S1_avg,'b')
        plot(1:length(intensityB2S1_avg),intensityB2S1_avg,'r')
        title('Average Intensity Side')
       
        subplot(1,3,3)
        hold on
        plot(1:length(intensityS2_avg),intensityS2_avg,'g')
        plot(1:length(intensityB1S2_avg),intensityB1S2_avg,'b')
        plot(1:length(intensityB2S2_avg),intensityB2S2_avg,'r')
        title('Average Intensity Side')
       
        subplot(1,3,1); 
    end
    hold on
    plot(1:length(intensity_avg),intensity_avg,'g')
    plot(1:length(intensity_tavg),intensity_tavg,'k')
    plot(1:length(intensityB1_avg),intensityB1_avg,'b')
    plot(1:length(intensityB2_avg),intensityB2_avg,'r')
    title('Average Intensity'); xlabel('Frame');ylabel('Intensity')
end

% Kymograph
if (nkymo > 0)
    figure
    map = colormap(jet(255));
    map = vertcat([0 0 0],map);
    imshow(uint8(kymo_avg.*255),map);
end

if (details == 1)
    figure
    subplot(1,2,1)
    histogram(Chist1(Chist1>0.1))
    hold on; histogram(Chist2(Chist2>0.1))
    title('Histogram C')
    
    subplot(1,2,2)
    histogram(ChistF1(ChistF1>0.1))
    hold on; histogram(ChistF2(ChistF2>0.1))
    title('Histogram CF')
    
    figure
    subplot(1,2,1)
    histogram(B1hist1(B1hist1>0.1))
    hold on; histogram(B1hist2(B1hist2>0.1))
    title('Histogram B1')
    
    subplot(1,2,2)
    histogram(B1histF1(B1histF1>0.1))
    hold on; histogram(B1histF2(B1histF2>0.1))
    title('Histogram B1F')
    
    figure
    subplot(1,2,1)
    histogram(B2hist1(B2hist1>0.1))
    hold on; histogram(B2hist2(B2hist2>0.1))
    title('Histogram B2')
    
    subplot(1,2,2)
    histogram(B2histF1(B2histF1>0.1))
    hold on; histogram(B2histF2(B2histF2>0.1))
    title('Histogram B2F')
end

