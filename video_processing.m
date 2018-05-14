function L = video_processing(pathf,fname,stp,smp,timestep,L,Cmin,Cmin_tmp,Cmax)

movie = [pathf '/' fname '_ratio.avi'];

V = VideoWriter(movie);
V.FrameRate = 1;
open(V);    

for count = 1:size(L,3)
    h = figure;
    figure(h);
    map = colormap(jet(255));
    map = vertcat([0 0 0],map);
    disp(['Video Processing:' num2str((count+stp-1))]);
    imshow(L(:,:,count),map); 
    txtstr = strcat('Time(s): ',num2str((count+stp-1)*timestep));
    text(10,10,txtstr,'color','white')
    hcb = colorbar;
    set(hcb,'FontSize',20)
    set(hcb,'YTick', [0 255])
    set(hcb,'YTickLabel', {num2str(Cmin_tmp),num2str(Cmax)})
    frame = getframe(gcf);
    set(gca, 'CLim', [0,255]);
    writeVideo(V,frame);
    close(h);
end

close(V);
disp(['Cmax:' num2str(Cmax)]);
disp(['Cmin:' num2str(Cmin*Cmax)]);

