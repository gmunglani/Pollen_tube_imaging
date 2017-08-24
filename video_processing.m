function [Cmin, Cmax, BN1, BN2, M] = video_processing(movie,stp,smp,BF1,BF2,posfront,framerate,time_step,Cmax,Cmin)

V = VideoWriter(movie);
V.FrameRate = framerate;
open(V);
    
for count = 1:smp-stp+1
    count
    
    BN1 = imcrop(BF1(:,:,count),posfront);
    BN2 = imcrop(BF2(:,:,count),posfront);
        
    type = find_orient(BN1);
    if (type == 1) BN1 = imrotate(BN1,-90); BN2 = imrotate(BN2,-90);
    elseif (type == 3) BN1 = imrotate(BN1,90); BN2 = imrotate(BN2,90);
    elseif (type == 4) BN1 = imrotate(BN1,180); BN2 = imrotate(BN2,180);
    end
        
    C = BN1./BN2;
    C(C==Inf) = 0;
    C(isnan(C)) = 0;
    M(:,:,count) = C;
end

if (Cmax == 0)
    Cmax = max(M(:))
end

K = M./Cmax;
K(isnan(K)) = 0;

if (Cmin == 0)
    %     Cmin = min(K(K>min(K)));
    Ktemp = bsxfun(@gt, K, min(K)); % bsxfun for Matlab versions older than 2016b
    Cmin = min(K(Ktemp));
end

% L = (K-Cmin)./(1-Cmin);
L =  bsxfun(@rdivide, bsxfun(@minus, K, Cmin),bsxfun(@minus, 1, Cmin));
L(L<0) = 0;
L = uint8(L.*255);

map = colormap(jet(255));
map = vertcat([0 0 0],map);

for k = 1:size(L,3)
    imshow(L(:,:,k),map);
    txtstr = strcat('Frame: ',num2str(k*time_step));
    text(10,10,txtstr,'color','white')
    hcb = colorbar;
    set(hcb,'YTick',[]);
    frame = getframe(gcf);
    set(gca, 'CLim', [0,255]);
    writeVideo(V,frame);
end
close(V);
Cmin
Cmax
