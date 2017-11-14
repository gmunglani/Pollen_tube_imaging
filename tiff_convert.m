clear all
close all

impath = '~/Documents/MATLAB/images/Ratio_tubes/TIFF/';
fname_YFP = 'YC_1_YFP.tif';
A1 = imread([impath fname_YFP],100); 
[tmp,pos] = imcrop(imagesc(A1));

for count = 1:100
     A1 = imread([impath fname_YFP],count); 
     AF1 = imcrop(A1,pos);
     imwrite(AF1,[impath 'tester2/YC1_YFP_test' num2str(count) '.tif'])
end
