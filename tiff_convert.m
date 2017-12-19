clear all
close all

impath = '~/Documents/Work/Images/Ratio_tubes/TIFF/';
fname_YFP = 'YC_1_YFP.tif';
A1 = imread([impath fname_YFP],1); 
%[tmp,pos] = imcrop(imagesc(A1));

for count = 1:100
     A1 = imread([impath fname_YFP],count); 
  %   AF1 = uint16((double(A1)./4095).*65535);
   %  AB1 = imcrop(A1,pos);
     imwrite(A1,[impath 'Poster_YC_vid' num2str(count) '.tif']);
end
