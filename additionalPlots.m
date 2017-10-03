moviesec = [movie '_sec'];

for count = stp:smp
    a = boundvid(:,:,count); count
    bound = a(any(a,2),:);
    out = zeros(size(C,1),size(C,2));
    out(sub2ind(size(out),bound(:,1),bound(:,2))) = 1;
    for j=1:length(bound)-1
        out = drawline(out,bound(j,1),bound(j,2),bound(j+1,1),bound(j+1,2),1); 
    end
%    Mset = M(:,:,count).*Fvid(:,:,count); Mset(Mset>1) = 1; 
%    Mset2 = M(:,:,count); Mset2(Mset2>1) = 1;
    Mvid(:,:,count) = (M(:,:,count).*Fvid(:,:,count)) + 100.*out;
end

video_processing(moviesec,stp,smp,BT1,BT2,framerate,timestep,Cmax,Cmin,Mvid);

out =[];
% movie = [movie '_line.avi'];
% V = VideoWriter(movie);
% V.FrameRate = framerate;
% open(V);

for count = smp:smp
%     h = figure;
%     figure(h);
    out = zeros(size(C,1),size(C,2));
    linectvid(:,:,count) = round(linectvid(:,:,count));
    out(sub2ind(size(out),linectvid(:,1,count),linectvid(:,2,count))) = 1;
    for j=1:100
        out = drawline(out,linectvid(j,1,count),linectvid(j,2,count),linectvid(j+1,1,count),linectvid(j+1,2,count),1);
    end
    for k=1:length(xy1f)-1
        out = drawline(out,xy1f(k,1),xy1f(k,2),xy2f(k,1),xy2f(k,2),1); 
    end
    out = imcomplement(out);

    DF = uint8((M(:,:,count).*out).*255);
    imshow(DF);
%     frame = getframe(gcf);
%     writeVideo(V,frame);
%     close(h);
end

%close(V);

% movieA1 = [movie '_YFP.avi'];
% V = VideoWriter(movie);
% V.FrameRate = framerate;
% open(V);
% for count = stp:smp
%     h = figure;
%     figure(h);
%     A1 = imread(fname_YFP,count); 
%     imshow(A1, [0 255])
%     frame = getframe(gcf);
%     writeVideo(V,frame);
%     close(h);
% end
% close(V);
% 
% A2 = imread(fname_CFP,count);