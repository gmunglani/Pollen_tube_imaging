clear all
close all

% Path to Mat file
path = '/Users/htv/Desktop/Background_Analysis_Results'; % Input folder path
fname = 'YC_3'; % File name 
stp = 64; % Start frame number
smp = 64; % End frame number

% Options for analysis
tip_plot = 0; % Video tip detection
video_intensity = 0; % Video intensity
frame_rate = 0.5; % Number of seconds per frame of input video
distributions = 0;  % Show histogram of results in the end
workspace = 0; % Save workspace

% Tip detection parameters
weight = 0.1; % Distance to eliminate branches (Higher means more reliance on the tip ellipse)

% ROI options
ROItype = 1; % No ROI = 0; Moving ROI = 1; Stationary ROI = 2
split = 1; % Split ROI along center line
circle = 0; % Circle ROI as fraction of diameter
starti = 0; % Rectangle ROI Start length / no pixelsize means percentage as a fraction of length of tube
stopi = 10; % Rectangle/Circle ROI Stop length / no pixelsize means percentage as a fraction of length of tube
pixelsize = 0.17; % Pixel to um conversion

% Kymo, movie and measurements options
Cmin = 2.5; % Min pixel value
Cmax = 4; % Max pixel value
nkymo = 3; % Number of pixels line width average for kymograph (odd number) (0 means no kymo)
diamcutoff = 0; % In pixels if pixelsize is not given

% Input ratio matrix
pathf = [path '/' fname]; % File path
M = h5read([pathf '/' fname '_back_proc_bleach.h5'],'/M'); % Ratio files
BT1 = h5read([pathf '/' fname '_back_proc_bleach.h5'],'/BT1'); % YFP
BT2 = h5read([pathf '/' fname '_back_proc_bleach.h5'],'/BT2'); % CFP

% Scaled plot of the growing tube with tip and ROI
if (tip_plot == 1)
    V = VideoWriter([pathf '/' fname '_growth.avi']);
    V.FrameRate = 1;
    open(V);
end

if (nkymo > 0 || video_intensity > 0) 
K = M(:,:,:)./Cmax;
    K(isnan(K)) = 0;
    Cmin_tmp = Cmin;
    Cmin = Cmin/Cmax;

    L = bsxfun(@rdivide, bsxfun(@minus, K, Cmin),bsxfun(@minus, 1, Cmin));
    L(L<0) = 0;
    L = uint8(L.*255);
end

% Make a movie and output min and max intensities of the whole stack
if (video_intensity)
    video_processing(pathf,fname,stp,smp,frame_rate,L(:,:,stp:smp),Cmin,Cmin_tmp,Cmax);
end

% Loop backwards over stack
for count = smp:-1:stp
    disp(['Image Analysis:' num2str(count)]);
    O = M(:,:,count);
    P = imbinarize(O,0.5);
    se = strel('disk',10);
    se2 = strel('disk',1);
    U = imopen(P,se);
    U = bwareaopen(U,100);
    U = bwareafilt(P,1);
    U = bwmorph(U,'clean');
    U = medfilt2(U);
    U = imclose(U,se);
    
    if (count == smp) Ub = logical(ones(size(P)));
    else Ub = U;
    end
        
    Urat = and(U,Ub);
    Ucount(count) = nnz(Urat)/nnz(U);
    if (Ucount(count) < 0.95 || count == smp) last_flag = 0;
    else last_flag = 1;
    end
    
    Umax = max(find(U(:,end)==1)); Umin = min(find(U(:,end)==1));
    U = imfill(drawline(U,Umin,size(U,2),Umax,size(U,2),1),'holes');

    % Removing branches from thinned image
    Q = bwmorph(U,'thin',Inf);    
    
    Qe = bwmorph(Q,'endpoints');
    [Qer,Qec] = find(Qe > 0);
    Qel = [Qer Qec];
        
    Qb = bwmorph(Q,'branchpoints');       
    [Qbr,Qbc] = find(Qb > 0);
    if (Qbr > 0)
        Qbf = [Qbr Qbc];
        [Q2, Qef, tmp] = branch_removal(Q,Qbf,Qel,0,1);
    else 
        Q2 = Q;
        [tmp,Qepos] = max(Qec);
        Qef = Qel;
        Qef(Qepos,:) = [];
    end
   
    % Finding the radius for ellipse fitting
    tols = 0; rad=1;
    while (tols == 0)
        sides = false; connect = false;
        try
            K = U(Qef(1)-rad:Qef(1)+rad,Qef(2)-rad:Qef(2)+rad);
        catch
            rad = 100;
            break;
        end
        Ke = [K(:,1)' K(end,2:end) K(end-1:-1:1,end)' K(1,end-1:-1:2)];
        Kd = diff(Ke);
        Kd(end+1) = Ke(1) - Ke(end);
        if (nnz(Kd) == 2) connect = true; end
        Ks = [sum(K(:,1)) sum(K(1,:)) sum(K(:,end)) sum(K(end,:))];
        if (nnz(Ks) <= 2)
            bou = bwboundaries(K); Kb = bou{1};
            Kbl = [find(Kb(:,1) == 1); find(Kb(:,2) == 1); find(Kb(:,1) == size(K,1)); find(Kb(:,2) == size(K,1))];
            if (size(Kbl,1) < size(K,1)) sides = true; end
        end
        if(connect == true && sides == true) tols = rad; end
        rad = rad + 1;
    end

    [boundb, tip_ellipse, tip_new, tip_check, diam, maxy, center, phin, axes, stats, edges] = locate_tip(U, tols, Qef);
    tip_ellipsepos = dsearchn(boundb,tip_ellipse);
    tip_ellipsef = boundb(tip_ellipsepos,:);
   
    % Skeletonizing and finding endpoints
    S = bwmorph(U,'skel',Inf);
    Se = bwmorph(S,'endpoints');
    [Ser,Sec] = find(Se > 0);
    Sel = [Ser Sec];
    
    % Skeltonizing and finding branchpoints
    Sb = bwmorph(S,'branchpoints');
    [Sbr,Sbc] = find(Sb > 0);
    Sbl = [Sbr Sbc];
    
    % Find branch point closest to thin edge
    if (last_flag == 0) Sbf = Sbl(dsearchn(Sbl,Qef),:);
    else [tmp, Sbmin] = min(pdist2(Sbl,tip_final_last) + pdist2(Sbl,Qef));
        Sbf = Sbl(Sbmin,:);
    end
    
    % Try and remove branches within some parameters
    close_dist = 0; 
    if (count == smp) diamo = diam; end
    if (pdist2(Qef,Sbf) > weight*diamo) close_dist = 1; end
    [S2,Sef,S2area] = branch_removal(S,Sbf,Sel,75,close_dist);
        
    % If more than 2 branches, further evaluation is needed
    if (size(Sef,1) > 1)
        % Voting to decide which branch to be chosen as closer to the tip
        [tmp, skel_ellipsepos] = min(pdist2(Sef,tip_ellipsef));
        if (last_flag == 1)
            [tmp, skel_lastpos] = min(pdist2(Sef,tip_final_last));
            if ((skel_lastpos+skel_ellipsepos+1) > 4) choice = 2; else choice = 1; end
        else 
            choice = skel_ellipsepos;       
        end
    
        tip_choice = [dsearchn(boundb,Sef(1,:));dsearchn(boundb,Sef(2,:))];
        if (min(tip_choice) == tip_choice(2)) S2area = 1/S2area; end
        tip_skelpos = tip_choice(choice);
        tip_skel = boundb(tip_skelpos,:);
   
        % Finding the middle of both branches if they exist
        cn = 0; tip_angle = [];
        for i = min(tip_choice):max(tip_choice)
            cn = cn+1;
            tip_angle(cn) = atan2((boundb(i,2) - Sbf(2)),(boundb(i,1) - Sbf(1)));
            if (pi - abs(max(tip_angle)) < abs(min(tip_angle)))
                if (tip_angle(cn) < 0) tip_angle(cn) = 2*pi + tip_angle(cn); end
            end
        end
        target_angle = (tip_angle(1) + tip_angle(end)*S2area)/(S2area+1);

        tip_anglediff = abs(tip_angle - target_angle);
        [tmp, tip_anglepos] = min(tip_anglediff);
        tip_midpos = tip_anglepos+min(tip_choice)-1;
        tip_mid = boundb(tip_midpos,:);
        
        test = [];
        if (tip_ellipsepos>min(tip_choice) && tip_ellipsepos<max(tip_choice))
            tip_final(count,:) = tip_ellipsef;
        else
            tip_ellipsedist = [pdist2(tip_ellipsef,tip_mid) pdist2(tip_ellipsef,tip_skel)];
            if (last_flag)            
                tip_finaldist = [pdist2(tip_final_last,tip_mid) pdist2(tip_final_last,tip_skel)];             
                [tmp, tip_finalpos] = min([(1-0.33)*tip_finaldist(1)+0.33*tip_ellipsedist(1) (1-0.33)*tip_finaldist(2)+0.33*tip_ellipsedist(2)]);
            else
                [tmp, tip_finalpos] = min(tip_ellipsedist);
            end
            if (tip_finalpos == 1) tip_final(count,:) = tip_mid; else tip_final(count,:) = tip_skel; end
            test = [test count];
        end
    else
        tip_skel = boundb(dsearchn(boundb,Sef(1,:)),:);
        if (last_flag) [tmp, tip_finaldistpos] = min([pdist2(tip_final_last,tip_ellipsef) pdist2(tip_final_last,tip_skel)]);   
        else tip_finaldistpos = 2;
        end
        if (tip_finaldistpos == 1) tip_final(count,:) = tip_ellipsef; else tip_final(count,:) = tip_skel; end
    end
        
    % Update tip_final for the next frame
    tip_final_last(:,:) = tip_final(count,:);
    
    % Find the curves along the sides of the tubes
    total1 = []; total2 = [];
    range1 = ceil(length(boundb)*0.5):length(boundb);
    dist1 = pdist2(boundb(range1,:),tip_final(count,:));
    postotal1 = find(dist1 > diamo*0.75)+range1(1)-1;
    if (~isempty(find(diff(postotal1(1:floor(length(postotal1)/2))>1))))
        postotal1(1:find(diff(postotal1(1:floor(length(postotal1)/2))>1))) = [];
    end
    total1(:,:) = boundb(postotal1,:);
    
    range2 = ceil(length(boundb)*0.5)-1:-1:1;
    dist2 = pdist2(boundb(range2,:),tip_final(count,:));
    postotal2 = range2(1)-find(dist2 > diamo*0.75)+1;
    if (~isempty(find(diff(postotal2(1:floor(length(postotal2)/2))>1))))
        postotal2(1:find(diff(postotal2(1:floor(length(postotal2)/2))>1))) = [];
    end
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
    
    % Find center line
    Q2line = drawline(Q2,tip_final(count,1),tip_final(count,2),Qef(1),Qef(2),1);  
    Q2bline1 = bwmorph(Q2line,'branchpoints');
    if (nnz(Q2bline1) > 0) 
        Q2line = imclose(Q2line,se2); 
        Q2line = bwmorph(Q2line,'thin'); 
    end
    
    Sbdist = pdist2(Sbl,[edges; round(mean(edges(:,1))),edges(2,2)]);
    [tmp, Sbpos] = min(Sbdist(:,3).*Sbdist(:,1)./Sbdist(:,2));
           
    Q2line = drawline(Q2line,round(mean(edges(:,1))),edges(2,2),Sbl(Sbpos,1),Sbl(Sbpos,2),1);
    fuse = 0;
    while(fuse == 0)
        Q2lineo = bwconncomp(Q2line);
        if (Q2lineo.NumObjects > 1)
            Q2line = imdilate(Q2line,se2);
        else
            fuse = 1;
            Q2line = bwmorph(Q2line,'thin',Inf);    
            break;
        end
    end
    
    Q2bline = bwmorph(Q2line,'branchpoints');
    [Qbrline,Qbcline] = find(Q2bline > 0);
    
    Q2eline = bwmorph(Q2line,'endpoints');
    [Qerline,Qecline] = find(Q2eline > 0);
    
    xct = []; yct = []; xc =[]; yc =[]; distct = []; distc = [];
    [Q3, tmp, tmp] = branch_removal(Q2line,[Qbrline Qbcline],[Qerline Qecline],0,2);
    Q3eline = bwmorph(Q3,'endpoints');
    [Q3erline,Q3ecline] = find(Q3eline > 0);
    [tmp,Q3pos] = min(Q3ecline);
    Q3geo = bwdistgeodesic(logical(Q3),Q3ecline(Q3pos),Q3erline(Q3pos),'quasi-euclidean');
    Q3geo(isnan(Q3geo)) = 0;
    [yct, xct, distct] = find(Q3geo);
    [distct geo_order] = sort(distct);
    yctk = yct(geo_order,:); xctk = xct(geo_order,:); 
    nline = 1:100; norder = floor(nline*max(distct)/100); nfinal = dsearchn(distct,norder');
    
    yct = yctk(nfinal); xct = xctk(nfinal); distct = distct(nfinal);
    xct = round(sgolayfilt(xct,3,15)); yct = round(sgolayfilt(yct,3,15)); 
    distc_t = pdist2(tip_final(count,:),Qef);
    [tmp, cut] = min(abs(distct - distc_t));
    xc = xct(cut:end); yc = yct(cut:end); distc = distct(cut:end); 
    
    % Calculate the gradient of the center line to get the normals
    dx = gradient(xc);
    dy = gradient(yc);
    
    % Finding the points where the normals hit the edge curves
    poscross1 = []; poscross2 = [];
    for n = 1:length(xc)
        nfitc = polyfit(vertcat(xc(n),(xc(n) - dy(n))),vertcat(yc(n),(yc(n) + dx(n))),1);
        if (n == 1) start_nfitc(:,:) = nfitc(:,:); end
        
        edge1 = total1(:,1) - nfitc(1).*total1(:,2) - nfitc(2);
        [tmp cross1] = min(abs(edge1));
        poscross1(n) = cross1;
        
        edge2 = total2(:,1) - nfitc(1).*total2(:,2) - nfitc(2);
        [tmp cross2] = min(abs(edge2));
        poscross2(n) = cross2;
    end
    
    % Ensure that all overlapping diameter lines are shifted backwards to
    % ensure continuity     
    [poscross1, poscross2, distcf] = line_continuity(poscross1,poscross2,1,distc);
    [poscross1, poscross2, distcf] = line_continuity(poscross1,poscross2,2,distcf);
     
    xy1 = []; xy2 = []; xy1 = total1(poscross1,:); xy2 = total2(poscross2,:); 
    if (length(xy1) > 20)
        xy1 = floor(sgolayfilt(xy1,3,15)); xy2 = floor(sgolayfilt(xy2,3,15));
    end
    xyout = vertcat(find(xy1(:,2) > size(U,2)), find(xy2(:,2) > size(U,2)));
    xy1(xyout,:) = []; xy2(xyout,:) = []; distcf(xyout) = [];
    
    % Cutoff the tip 
    [tmp, distpos, tmp] = intersect(distc,distcf);
    distctf = [distct(1:cut-1); distc(distpos)]; xctf = [xct(1:cut-1); xc(distpos)]; yctf = [yct(1:cut-1); yc(distpos)]; 
    linectf = [yctf xctf];

    if (ROItype > 0)
        Esize = size(U);
        % Find ROI from centerline distance using percentages or distance
        if (pixelsize == 0)
            percent = (100*distctf)./(distctf(end));
            start_length = abs(percent - starti); [tmp startpos] = min(start_length);
            stop_length = abs(percent - stopi); [tmp stoppos] = min(stop_length);
            distc_t = (100*distc_t)/max(distctf);
        else
            start_length = abs(distctf*pixelsize - starti); [tmp startpos] = min(start_length);
            stop_length = abs(distctf*pixelsize - stopi); [tmp stoppos] = min(stop_length);
            distc_t = distc_t*pixelsize;
        end
        
        % Project ROI length onto the side curves
        [startc1,stopc1] = closest_bound(total1,xctf,yctf,max(startpos,1),max(stoppos,1));
        [startc2,stopc2] = closest_bound(total2,xctf,yctf,max(startpos,1),max(stoppos,1));
        
        % Create masks for rectangles and circles, and include whether they are
        % normal, split or stationary
        if (ROItype ~= 2 | count == smp)
            if (circle == 0)
                roi = vertcat(total1(startc1:stopc1,:), total2(stopc2:-1:startc2,:));
                if (starti < distc_t) roi = vertcat(boundb(postotal2(1):postotal1(2),:),roi); end
                F = poly2mask(roi(:,2),roi(:,1),Esize(1),Esize(2));
            else
                mask = zeros(Esize(1),Esize(2));
                roi = [linectf(stoppos,1) linectf(stoppos,2)];
                mask(roi(1),roi(2)) = 1;
                F = bwdist(mask) >= 0.5*circle.*diamo;
                F = imcomplement(F);
            end
            
            if (split == 1)
                if (circle > 0)
                    stoppos = length(linectf); stopc1 = length(total1); stopc2 = length(total2);
                end
                roi1 = vertcat(total1(startc1:stopc1,:), linectf(stoppos:-1:startpos,:));
                roi2 = vertcat(total2(startc2:stopc2,:), linectf(stoppos:-1:startpos,:));
                if (starti < distc_t)
                    roi1 = vertcat(boundb(range2(1):postotal1(2),:),roi1,boundb(range2(1),:));
                    roi2 = vertcat(boundb(range2(1):-1:postotal2(1),:),roi2,boundb(range2(1),:));
                end
                F1 = F.*poly2mask(roi1(:,2),roi1(:,1),Esize(1),Esize(2));
                F2 = F.*poly2mask(roi2(:,2),roi2(:,1),Esize(1),Esize(2));
            end
        end
    
    
        % Calculate average intensities and pixel numbers
        Fpixelnum(count) = nnz(O.*F);
        intensityM(count) = sum(O(:))/nnz(O);
        intensityM_F(count) = sum(sum(O.*F))/Fpixelnum(count);
        intensityB1_F(count) = sum(sum(BT1(:,:,count).*F))/Fpixelnum(count);
        intensityB2_F(count) = sum(sum(BT2(:,:,count).*F))/Fpixelnum(count);
        if (split)
            F1pixelnum(count) = nnz(O.*F1);
            F2pixelnum(count) = nnz(O.*F2);
            intensityM_F1(count) = sum(sum(O.*F1))/F1pixelnum(count);
            intensityB1_F1(count) = sum(sum(BT1(:,:,count).*F1))/F1pixelnum(count);
            intensityB2_F1(count) = sum(sum(BT2(:,:,count).*F1))/F1pixelnum(count);
            intensityM_F2(count) = sum(sum(O.*F2))/F2pixelnum(count);
            intensityB1_F2(count) = sum(sum(BT1(:,:,count).*F2))/F2pixelnum(count);
            intensityB2_F2(count) = sum(sum(BT2(:,:,count).*F2))/F2pixelnum(count);
        end
        
        % Histogram of first and last frame
        if (distributions)
            d = 1;
            if (count == stp || count == smp)
                Msize = [numel(O),1]; BT1size = [numel(BT1(:,:,count)),1]; BT2size = [numel(BT2(:,:,count)),1];
                Mhist(:,d) = reshape(O,Msize);
                B1hist(:,d) = reshape(BT1(:,:,count),BT1size);
                B2hist(:,d) = reshape(BT2(:,:,count),BT2size);
                
                MhistF(:,d) = reshape(O.*F,Msize);
                B1histF(:,d) = reshape(BT1(:,:,count).*F,BT1size);
                B2histF(:,d) = reshape(BT2(:,:,count).*F,BT2size);
            end
            d = d+1;
        end
    end
    
    % Cut off the tip part of the diameter calculation if necessary
    if (pixelsize > 0) cutoffp = dsearchn(distcf',diamcutoff/pixelsize);
    else cutoffp = dsearchn(distcf',diamcutoff);
    end
    if (cutoffp > 1) xy1(cutoffp-1,:) = []; xy2(cutoffp-1,:) = []; end

    % Diameter of tube
    diamf = diag(pdist2(xy1,xy2));
    diamf_avg(count) = sum(diamf)/length(diamf);
    
    % Kymograph
    if (nkymo > 0)
        if (count == smp) npoints = round(distct(end)*1.1); end
        
        % Average number of points across the tube width in kymo based on orientation
        linecte = []; linecte(:,:,1) = [yctk, xctk];
        for a = 2:nkymo
            if (mod(a,2) == 0)
                ind = floor(a*0.5);
            else
                ind = -floor(a*0.5);
            end
            if (start_nfitc(1) < 0)
                if (mod(a,2) == 0) linecte(:,:,a) = [yctk+ind, xctk-ind];
                else linecte(:,:,a) = [yctk-ind, xctk+ind];
                end
            else
                if (mod(a,2) == 0) linecte(:,:,a) = [yctk+ind, xctk+ind];
                else linecte(:,:,a) = [yctk-ind, xctk-ind];
                end
            end
        end
        kymo = [];
        for a = 1:nkymo
            kymo(:,a) = improfile(imgaussfilt(L(:,:,count),1.5), linecte(:,2,a), linecte(:,1,a), double(round(distct(end))));
        end
        kymo(isnan(kymo)) = 0;
        kymo_avg(:,count-stp+1) = vertcat(zeros((5 + npoints - round(distct(end))),1), mean(kymo,2));
    end 

    % Tip plot
    Splot = zeros(size(Q2));
    Splot(tip_final(count,1)-3:tip_final(count,1)+3,tip_final(count,2)-3:tip_final(count,2)+3) = 1;
 %   Splot(tip_ellipsef(1)-1:tip_ellipsef(1)+1,tip_ellipsef(2)-1:tip_ellipsef(2)+1) = 2;
 %   if (size(Sef,1) > 1) Splot(tip_mid(1)-3:tip_mid(1)+3,tip_mid(2)-3:tip_mid(2)+3) = 3; end
 %   Splot(tip_skel(1)-1:tip_skel(1)+1,(2)-1:tip_skel(2)+1) = 4;
    
    Cplot = zeros(size(Q2)); Cplot(sub2ind([size(Cplot,1) size(Cplot,2)],yctk,xctk)) = 2.*ones(size(xctk));
    for j = 1:length(xy1) Cplot = drawline(Cplot,xy1(j,1),xy1(j,2),xy2(j,1),xy2(j,2),1); end
    Cplot(:,size(U,2)+1:end) = [];
    
    % Plot two images
    if (tip_plot) h = figure('visible', 'off');
    else h = figure;
    end
    
    subplot(1,2,1)
    image1 = U+S2*2+Splot*4;
    imagesc(image1)
    
    subplot(1,2,2)
    image2 = U+Cplot.*2;
    if (ROItype > 0) image2 = image2 +F*4; end
    imagesc(image2);
    
    
    if (tip_plot)
        txtstr = strcat('Frame: ',num2str((count)));
        text(10,10,txtstr,'color','white')
        txtstr = strcat('Time(s): ',num2str((count*frame_rate)));
        text(50,10,txtstr,'color','white')
        frame = getframe(gcf);
        writeVideo(V,frame);
        close(h);
    end
end

if (tip_plot == 1) close(V); end

% Final tip movement/diameter/pixel number on a per frame basis
figure
subplot(3,1,1)
plot(tip_final(stp:smp,2),tip_final(stp:smp,1),'b')
axis([min(tip_final(stp:smp,2))-5  max(tip_final(stp:smp,2))+5 min(tip_final(stp:smp,1))-5 max(tip_final(stp:smp,1))+5]);
title('Tip Final Position', 'FontSize',16);

subplot(3,1,2)
plot(stp:smp,diamf_avg(stp:smp),'b')
xlabel('Frame', 'FontSize',12);
title('Average Diameter','FontSize',16)
axis([stp-1 smp+1 0.5*max(diamf_avg) 1.25*max(diamf_avg)])

subplot(3,1,3)
plot(stp:smp,intensityM(stp:smp),'k') % Ratio only
hold on
plot(stp:smp,intensityM_F(stp:smp),'r') % Ratio and full ROI
if (split)
    plot(stp:smp,intensityM_F1(stp:smp),'b*') % Ratio and split ROI 1
    plot(stp:smp,intensityM_F2(stp:smp),'g*') % Ratio and split ROI 2
end
xlabel('Frame', 'FontSize',12);
title('Intensity Ratio', 'FontSize',16)
axis([stp-1 smp+1 1 max(intensityM)*1.25]);

% Kymograph
if (nkymo > 0)
    kymo_avg(find(kymo_avg<0)) = 0;
    figure
    map = colormap(jet(255));
    map = vertcat([0 0 0],map);
    imshow(uint8(kymo_avg.*255/max(kymo_avg(:))),map);
end

% Total intensity plots
if (ROItype > 0)
    figure
    if (split) 
        F1ratio = intensityB1_F1(stp:smp)./intensityB2_F1(stp:smp);
        subplot(1,3,2)
        hold on
        plot(stp:smp,F1ratio,'b')
        axis([stp-1 smp+1 0.8 max(F1ratio(:))*1.25]);
        title('Intensity F1 (split ROI 1)'); xlabel('Frame');
       
        F2ratio = intensityB1_F2(stp:smp)./intensityB2_F2(stp:smp);
        subplot(1,3,3)
        hold on
        plot(stp:smp,F2ratio,'b')
        axis([stp-1 smp+1 0.8 max(F1ratio(:))*1.25]);
        title('Intensity F2 (split ROI 2)'); xlabel('Frame');
       
        subplot(1,3,1); 
        plot(stp:smp,F2ratio./F1ratio,'b')
        axis([stp-1 smp+1 0.5 2]);
        title('Intensity ratio between split ROIs'); xlabel('Frame');
    else
        Fratio = intensityB1_F(stp:smp)./intensityB2_F(stp:smp);
        hold on
        plot(stp:smp,Fratio,'b');
        axis([stp-1 smp+1 0.8 max(Fratio(:))*1.25]);
        title('Intensity F'); xlabel('Frame');
    end
end

% Distributions of intensity on the first and last frames
if (distributions == 1)    
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

if (workspace) save([pathf '/' fname '_result.mat']); end