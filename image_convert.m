clear all
close all

fname = 'Data31_YFP';
fname_full = strcat('~/Documents/MATLAB/images/',fname,'.tif');
info = imfinfo(fname_full);
mImage=info(1).Width;
nImage=info(1).Height;
num_images = numel(info);

%modify=zeros(nImage,mImage,num_images,'uint16');
modify=zeros(nImage,mImage,num_images,'uint8');
for count=1:360
    inp = imread(fname_full,count);
   % modify(:,:,count) = uint16(double(inp)*65535/4096);
    modify(:,:,count) = uint8(double(inp)*255/4096);
  %  imwrite(modify(:,:,count), strcat('./images/',fname,'_images8/',fname,sprintf('%03d',count),'.tif'))
    imwrite(modify(:,:,count), strcat('~/Documents/MATLAB/images/Pollentube_images/image',sprintf('%03d',count),'.tif'))
end