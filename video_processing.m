function [L,Cmin, Cmax] = video_processing(movie,stp,smp,fname,fname2,pos,posback,thresh1,time_step,decay,gauss,framerate,Cmax,Cmin)

V = VideoWriter(movie);
V.FrameRate = framerate;
open(V);
    
for count = stp:smp
    count 
    A1 = imread(fname, count, 'Info', imfinfo(fname)); A1c = imcrop(A1,pos); A1c = A1c ./ exp(decay*time_step*count);    
    A2 = imread(fname2, count, 'Info', imfinfo(fname2)); A2c = imcrop(A2,pos);

    B1 = imgaussfilt(mat2gray(A1),gauss); B1c = imgaussfilt(mat2gray(A1c),gauss); 
    B1b = imcrop(B1,posback); B1b(B1b>thresh1) = 0; back1 = sum(sum(B1b))/length(find(B1b>0));
    
    B2 = imgaussfilt(mat2gray(A2),gauss); B2c = imgaussfilt(mat2gray(A2c),gauss); 
    B2b = imcrop(B2,posback); B2b(B2b>thresh1) = 0; back2 = sum(sum(B2b))/length(find(B2b>0));
    
    B1c(B1c<thresh1) = 0; B1n = B1c - back1; B1n(B1n<0) = 0;
    B2c(B2c<thresh1) = 0; B2n = B2c - back2; B2n(B2n<0) = 0;
    
    C = B1n./B2n;
    C(C==Inf) = 0;
    C(isnan(C)) = 0;
    M(:,:,count) = C;
end

if (Cmax == 0)
    Cmax = max(M(:));
end

K = M./Cmax;
K(isnan(K)) = 0;

if (Cmin == 0)    
    Cmin = min(K(K>min(K)));
end

L = (K-Cmin)./(1-Cmin);
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
