function [Cmin, Cmax] = video_processing(movie,stp,smp,BT1,BT2,framerate,timestep,Cmax,Cmin,M)

V = VideoWriter(movie);
V.FrameRate = framerate;
open(V);    

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
    txtstr = strcat('Frame: ',num2str(k*timestep));
    text(10,10,txtstr,'color','white')
    hcb = colorbar;
    set(hcb,'YTick',[]);
    frame = getframe(gcf);
    set(gca, 'CLim', [0,255]);
    writeVideo(V,frame);
end
close(V);
%Cmin
%Cmax
