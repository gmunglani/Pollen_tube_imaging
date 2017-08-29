function video_processing(movie,framerange,BT1,BT2,framerate,timestep,Cmax,Cmin,M)

movie = [movie '.avi'];
V = VideoWriter(movie);
V.FrameRate = framerate;
open(V);    

K = M(:,:,1:framerange)./Cmax;
K(isnan(K)) = 0;
Cmin = Cmin/Cmax;

% L = (K-Cmin)./(1-Cmin);
L = bsxfun(@rdivide, bsxfun(@minus, K, Cmin),bsxfun(@minus, 1, Cmin));
L(L<0) = 0;
L = uint8(L.*255);

map = colormap(jet(255));
map = vertcat([0 0 0],map);

figure
for count = 1:size(L,3)
    imshow(L(:,:,count),map);
    txtstr = strcat('Frame: ',num2str(count*timestep));
    text(10,10,txtstr,'color','white')
    hcb = colorbar;
    set(hcb,'YTick',[]);
    frame = getframe(gcf);
    set(gca, 'CLim', [0,255]);
    writeVideo(V,frame);
end
close(V);
disp(['Cmax:' num2str(Cmax)]);
disp(['Cmin:' num2str(Cmin*Cmax)]);

