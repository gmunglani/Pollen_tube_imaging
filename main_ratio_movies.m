clear all
close all

% Path to h5 file
path = '/home/gm/Documents/Scripts/MATLAB/Tip_results'; % Input folder path
fname = 'YC11'; % File name 
stp = 1; % Start frame number
smp = 200; % End frame number
specific = []; % Frames to change

% Bleach options
bleachYFP = 1:1; % Bleaching range YFP (Greater than length 1 commences bleaching)
bleachCFP = 1:1; % Bleaching range CFP (Greater than length 1 commences bleaching)

% Other Options
register = 1; % Register image
union = 1; % Take the union of the two image masks
mask_plot = 1; % Plot the mask and overlap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frame range
pathf = [path '/' fname '/' fname]; flag = 0;
if (exist([pathf '_back_proc.h5'],'file') == 2)
    if ((length(bleachYFP)+length(bleachCFP) > 2) || ~isempty(specific))
        BT1 = h5read([pathf '_back_proc.h5'],'/BT1');
        BT2 = h5read([pathf '_back_proc.h5'],'/BT2');
        M = h5read([pathf '_back_proc.h5'],'/M');
        intensity1 = h5read([pathf '_back_proc.h5'],'/intensity1');
        intensity2 = h5read([pathf '_back_proc.h5'],'/intensity2');
        Bsum = h5read([pathf '_back_proc.h5'],'/Bsum');
        
        posfront = h5readatt([pathf '_back_proc.h5'],'/M','Crop');
        Bedge = h5readatt([pathf '_back_proc.h5'],'/M','Bedge');
        mask_plot = 0;
    end

    % Bleach correction
    if (length(bleachYFP)+length(bleachCFP) > 2)
        disp('Bleach correction'); flag = 1;
        
        if(length(bleachYFP)>1)
            fit1 = fit(bleachYFP',intensity1(bleachYFP)','exp1','StartPoint',[intensity1(1),0.005]);
            for i = bleachYFP(1):smp
                BT1(:,:,i) = BT1(:,:,i).*intensity1(1)/fit1(i);
            end
        end
        
        if(length(bleachCFP)>1)
            fit2 = fit(bleachCFP',intensity2(bleachCFP)','exp1','StartPoint',[intensity2(1),0.005]);
            for i = bleachCFP(1):smp
                BT2(:,:,i) = BT2(:,:,i).*intensity2(1)/fit2(i);
            end
        end
    end
end
  
if (length(bleachYFP)+length(bleachCFP) == 2)
    if (mask_plot == 1)
        V2 = VideoWriter([pathf '_mask.avi']);
        V2.FrameRate = 1;
        open(V2);
    end
    
    YFP = h5read([pathf '_back.h5'],'/YFP');
    CFP = h5read([pathf '_back.h5'],'/CFP');
    
    if (isempty(specific))
        % Crop region on the last frame
        AC = YFP(:,:,smp);
        BC = mat2gray(AC);
        [tmp,posfront] = imcrop(BC');
        range = stp:smp;
        Bedge = zeros(1,smp); Bsum = zeros(1,smp);
    else
        range = specific;
    end
    
    [optimizer, metric] = imregconfig('multimodal');
    
    for count = range
        % Read image and add bleach correction
        disp(['Pre Processing:' num2str(count)]);
        
        A1 = YFP(:,:,count)';
        A2 = CFP(:,:,count)';
        
        B1 = imcrop(A1,posfront);
        B2 = imcrop(A2,posfront);

        if (register == 1) B2 = imregister(B2,B1,'translation',optimizer,metric); end
        
        Bu1 = uint8(single(B1).*255/4095);
        Bu2 = uint8(single(B2).*255/4095);
        
        level1 = graythresh(Bu1); Bl1 = imbinarize(Bu1,level1); 
        level2 = graythresh(Bu2); Bl2 = imbinarize(Bu2,level2);
        
        B1 = B1.*uint16(Bl1);
        B2 = B2.*uint16(Bl2);
        
        % Orient image
        if (count==range(1)) type = find_orient(B1); end
        if (type == 1) B1 = imrotate(B1,-90); B2 = imrotate(B2,-90);
        elseif (type == 3) B1 = imrotate(B1,90); B2 = imrotate(B2,90);
        elseif (type == 4) B1 = imrotate(B1,180); B2 = imrotate(B2,180);
        end
        
        % Union of images to ensure perfect overlap
        B = B1.*B2; B(B>0) = 1;
        
        % Analysis
        Bsum(count,1) = sum(B(:));
        Bsum(count,2) = nnz(B1);
        Bsum(count,3) = nnz(B2);
        
        if (isempty(specific))
            if (count == range(1))
                BT1=zeros(size(B1,1),size(B1,2),range(end));
                BT2=zeros(size(B2,1),size(B2,2),range(end));
            end
            while(sum(B(:,end-Bedge(count))) == 0)
                Bedge(count) = Bedge(count) + 1;
            end
        else
            B = B(:,1:end-Bedge);
            B1 = B1(:,1:end-Bedge);
            B2 = B2(:,1:end-Bedge);
        end
        
        if (union)
            BR1 = B1.*B;
            BR2 = B2.*B;
        else
            BR1 = B1;
            BR2 = B2;
        end
        
        BT1(:,:,count) = BR1; intensity1(count) = double(median(nonzeros(BR1(:))));
        BT2(:,:,count) = BR2; intensity2(count) = double(median(nonzeros(BR2(:))));
        
        
        if (mask_plot)
            h2 = figure('visible', 'off');
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
            plot(stp:count,intensity1(stp:count)/intensity1(stp),'b');
            hold on
            plot(stp:count,intensity2(stp:count)/intensity2(stp),'r');
            title(['Median Pixel Intensity:' num2str(count)]);
            xlabel('Frame')
            axis tight
            
            frame = getframe(gcf);
            writeVideo(V2,frame);
            close(h2);
        end
    end
    if (mask_plot == 1) close(V2); end
    if (isempty(specific))
        BT1 = BT1(:,1:end-max(Bedge),:);
        BT2 = BT2(:,1:end-max(Bedge),:);
    end
end

% Create ratio image and output values
M = BT1./BT2;
M(M==Inf) = 0;
M(isnan(M)) = 0;

% Write h5 file and delete old files
if (flag) name = [pathf '_back_proc_bleach.h5'];
else name = [pathf '_back_proc.h5'];
end

delete(name);
h5create(name,'/M',size(M));
h5create(name,'/BT1',size(BT1));
h5create(name,'/BT2',size(BT2));
h5create(name,'/intensity1',size(intensity1));
h5create(name,'/intensity2',size(intensity2));
h5create(name,'/Bsum',size(Bsum));

h5write(name,'/M',uint16(M));
h5writeatt(name,'/M','Crop',posfront);
h5writeatt(name,'/M','Bedge',max(Bedge));

h5write(name,'/BT1',uint16(BT1));
h5write(name,'/BT2',uint16(BT2));
h5write(name,'/intensity1',intensity1);
h5write(name,'/intensity2',intensity2);
h5write(name,'/Bsum',Bsum);

% Plotting signal results and decay
h = figure;
subplot(2,2,1)
plot(stp:smp,Bsum(stp:smp,1),'b');
hold on
plot(stp:smp,Bsum(stp:smp,2),'r');
plot(stp:smp,Bsum(stp:smp,3),'g');
axis([stp-1 smp+1 0.5*max(Bsum(:,1)) max(Bsum(:,1))])
title('Total Number of Pixels');

% Scale each channel by maximum intensity for viewing
BT1max = max(BT1(:));
BT1 = BT1./BT1max;
BT2 = BT2./BT1max;
    
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
axis([stp-1 smp+1 0 max(Mmax)])
set(gca,'YMinorTick','on')
title('Percentiles of Ratio image');

subplot(2,2,3)
hold on
plot(stp:smp,B1max(stp:smp),'b*')
plot(stp:smp,B1max_prc(stp:smp),'bd')
plot(stp:smp,B1min(stp:smp),'r*')
plot(stp:smp,B1min_prc(stp:smp),'rd')
grid on
axis([stp-1 smp+1 0 1])
set(gca,'YMinorTick','on');
title('Percentiles of YFP image');

subplot(2,2,4)
hold on
plot(stp:smp,B2max(stp:smp),'b*')
plot(stp:smp,B2max_prc(stp:smp),'bd')
plot(stp:smp,B2min(stp:smp),'r*')
plot(stp:smp,B2min_prc(stp:smp),'rd')
grid on
axis([stp-1 smp+1 0 1])
set(gca,'YMinorTick','on');
title('Percentiles of CFP image');
hold off

if (flag) savefig(h,[pathf '_back_proc_bleach.fig']);
else savefig(h,[pathf '_back_proc.fig']);
end

