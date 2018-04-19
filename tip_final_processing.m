clear all
close all

% Path to Mat file
path = '/Users/htv/Desktop/1403'; % Input folder path
fname = 'YC_2'; % File name 
stp = 1; % Start frame number
smp = 1802; % End frame number

% Bleach options
bleach1 = 1:1; % Bleaching range YFP
bleach2 = 1:1; % Bleaching range CFP

% Other Options
register = 1; % Register image
union = 1; % Take the union of the two image masks
mask_plot = 1; % Plot the mask and overlap
analysis = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frame range
pathf = [path '/' fname];
if  exist([pathf '/' fname '_proc.mat'],'file') == 2    
    load([pathf '/' fname '_proc.mat']);
    analysis = 0;
    disp('Analysis file exists');
else
    fname_YFP = h5read([pathf '/' fname '_YFP.h5'],'/values');
    fname_CFP = h5read([pathf '/' fname '_CFP.h5'],'/values');
 
    % Crop region on the last frame
    AC = fname_YFP(:,:,smp);
    BC = mat2gray(AC);
    [tmp,posfront] = imcrop(BC');
    [optimizer, metric] = imregconfig('multimodal');

    Bedge = zeros(1,smp); Bsum = zeros(1,smp);
    
    if (mask_plot == 1)
        V2 = VideoWriter([pathf '/' fname '_mask.avi']);
        V2.FrameRate = 1;
        open(V2);
    end
    
    for count = stp:smp
        % Read image and add bleach correction
        disp(['Pre Processing:' num2str(count)]);
        
        A1 = fname_YFP(:,:,count); A1n = A1';
        A2 = fname_CFP(:,:,count); A2n = A2';
        
        B1 = imcrop(A1n,posfront);
        B2 = imcrop(A2n,posfront);
        
        if (register == 1) B2 = imregister(B2,B1,'translation',optimizer,metric); end
        
        Bu1 = uint8(single(B1).*255/4095);
        Bu2 = uint8(single(B2).*255/4095);
        
        level1 = graythresh(Bu1); Bl1 = imbinarize(Bu1,level1); 
        level2 = graythresh(Bu2); Bl2 = imbinarize(Bu2,level2);
        
        B1 = B1.*uint16(Bl1);
        B2 = B2.*uint16(Bl2);
        
        % Orient image
        if (count==stp) type = find_orient(B1); end
        if (type == 1) B1 = imrotate(B1,-90); B2 = imrotate(B2,-90);
        elseif (type == 3) B1 = imrotate(B1,90); B2 = imrotate(B2,90);
        elseif (type == 4) B1 = imrotate(B1,180); B2 = imrotate(B2,180);
        end
        
        if (count == stp)
            BT1=zeros(size(B1,1),size(B1,2),smp);
            BT2=zeros(size(B2,1),size(B2,2),smp);
        end
         
        % Union of images to ensure perfect overlap
        B = B1.*B2; B(B>0) = 1;
        while(sum(B(:,end-Bedge(count))) == 0)
            Bedge(count) = Bedge(count) + 1;
        end

        if (union == 1)
            BT1(:,:,count) = B1.*B;
            BT2(:,:,count) = B2.*B;
        else
            BT1(:,:,count) = B1;
            BT2(:,:,count) = B2;
        end
        
        % Analysis
        Bsum(count) = sum(B(:));
        B1sum(count) = nnz(B1);
        B2sum(count) = nnz(B2);
        
        Br1 = reshape(BT1(:,:,count),[numel(BT1(:,:,count)),1]);
        Br2 = reshape(BT2(:,:,count),[numel(BT2(:,:,count)),1]);
        
        intensity1(count) = median(Br1);
        intensity2(count) = median(Br2);
             
        h2 = figure('Visible', 'off');
        subplot(2,2,1)
        hold on
        imshow(B1, [0 max(B1(:))]);
        subplot(2,2,2)
        hold on
        imshow(B2, [0 max(B2(:))]);
        subplot(2,2,3)
        hold on
        imshowpair(B1,B2);
        subplot(2,2,4)
        axis([stp smp 0 4096])
        plot(stp:count,intensity1(stp:count),'b');
        hold on
        plot(stp:count,intensity2(stp:count),'r');
        title(['Median Pixel Intensity:' num2str(count)]);
        xlabel('Frame')
        axis tight
        
        frame = getframe(gcf);
        writeVideo(V2,frame);
        close(h2);   
    end
    if (mask_plot == 1) close(V2); end

    if(length(bleach1)>1) 
        fit1 = fit(bleach1',intensity1(bleach1)','exp1','StartPoint',[intensity1(1),0.005]); 
        for i = bleach1(1):smp
            BT1(:,:,i) = BT1(:,:,i)*intensity1(1)/fit1(i);
        end
    end
    if(length(bleach2)>1) 
        fit2 = fit(bleach2',intensity2(bleach2)','exp1','StartPoint',[intensity2(1),0.005]); 
        for i = bleach2(1):smp
            BT2(:,:,i) = BT2(:,:,i)*intensity2(1)/fit2(i);
        end
    end
    
    BT1 = BT1(:,1:end-max(Bedge),:);
    BT2 = BT2(:,1:end-max(Bedge),:);
    
    BT1max = max(BT1(:))
    BT1 = BT1./BT1max;
    BT2 = BT2./BT1max; 
        
    % Create ratio image and output values
    M = BT1./BT2; 
    M(M==Inf) = 0;
    M(isnan(M)) = 0;
end

if (analysis == 1)
    h = figure;
    figure(h);
    subplot(2,2,1)
    plot(stp:smp,Bsum(stp:smp),'b');
    hold on
    plot(stp:smp,B1sum(stp:smp),'r');
    plot(stp:smp,B2sum(stp:smp),'g');
    axis([stp smp 0.5*max(Bsum) max(Bsum)])
    title('Total Number of Pixels');

    [Mmin Mmax Mmin_prc Mmax_prc] = channel_analysis(M,smp);
    [B1min B1max B1min_prc B1max_prc] = channel_analysis(BT1,smp);
    [B2min B2max B2min_prc B2max_prc] = channel_analysis(BT2,smp);
    
    subplot(2,2,2)
    hold on
    plot(stp:smp,Mmax(stp:smp),'b*')
    plot(stp:smp,Mmax_prc(stp:smp),'bd')
    plot(stp:smp,Mmin(stp:smp),'r*')
    plot(stp:smp,Mmin_prc(stp:smp),'rd')
    grid on
    axis([stp smp 0 max(Mmax)])
    set(gca,'YMinorTick','on')
    title('Percentiles of Ratio image');

    subplot(2,2,3)
    hold on
    plot(stp:smp,B1max(stp:smp),'b*')
    plot(stp:smp,B1max_prc(stp:smp),'bd')
    plot(stp:smp,B1min(stp:smp),'r*')
    plot(stp:smp,B1min_prc(stp:smp),'rd')
    grid on      
    axis([stp smp 0 1])
    set(gca,'YMinorTick','on');
    title('Percentiles of YFP image');
    
    subplot(2,2,4)
    hold on
    plot(stp:smp,B2max(stp:smp),'b*')
    plot(stp:smp,B2max_prc(stp:smp),'bd')
    plot(stp:smp,B2min(stp:smp),'r*')
    plot(stp:smp,B2min_prc(stp:smp),'rd')
    grid on
    axis([stp smp 0 1])
    set(gca,'YMinorTick','on');
    title('Percentiles of CFP image');
    hold off
    save([pathf '/' fname '_proc.mat'],'BT1','BT2','M');
    savefig(h,[pathf '/' fname '_proc.fig']);
    error('Check pixel counts and min/max intensity');
end
