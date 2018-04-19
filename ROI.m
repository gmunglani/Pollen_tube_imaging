
% Find ROI from centerline distance using percentages or distance
if (pixelsize == 0)
    percent = (distct)./(distct(end));
    start_length = abs(percent*100 - starti); [tmp startpos] = min(start_length)
    stop_length = abs(percent*100 - stopi); [tmp stoppos] = min(stop_length)
else
    start_length = abs(distct*pixelsize - starti); [tmp startpos] = min(start_length);
    stop_length = abs(distct*pixelsize - stopi); [tmp stoppos] = min(stop_length);
end

% Create masks for rectangles and circles, and include whether they are
% normal, split or stationary
if (ROItype ~= 2 | count == stp)
    if (circle == 0)
        % Project ROI length onto the side curves
        [startc1,stopc1] = closest_bound(total1,xy1(:,2),xy1(:,1),startpos,stoppos);
        [startc2,stopc2] = closest_bound(total2,xy2(:,2),xy2(:,1),startpos,stoppos);

        roi = vertcat(total1(startc1:stopc1,:), total2(stopc2:-1:startc2,:));
        if (starti < distc_t*pixelsize) roi = vertcat(boundb(postotal2(1):postotal1(2),:),roi); end
        F = poly2mask(roi(:,2),roi(:,1),Esize(1),Esize(2));
    else
        mask = zeros(Esize(1),Esize(2));
%         if (stopi < distc_t*pixelsize)
%             npts = round(distc_t/(distc(min(find(distc>1)))*pixelsize));
%             xcir = linspace(tip_final(:,2),round(linec(1,2)),npts);
%             ycir = linspace(tip_final(:,1),round(linec(1,1)),npts);
%             distcir = diag(pdist2([tip_final(:,1)*ones(npts,1) tip_final(:,2)*ones(npts,1)], [ycir; xcir]'));
%             
%             if (pixelsize > 0)
%                 stop_length = abs(distcir.*pixelsize - stopi);
%             else
%                 stop_length = abs(distcir - stopi);
%             end
%             [stoppos stoppos] = min(stop_length);
%             roi = [ycir(stoppos) xcir(stoppos)];
%             mask(round(ycir(stoppos)),round(xcir(stoppos))) = 1;
%         else
            roi = [round(linec(stoppos,1)) round(linec(stoppos,2))];
            mask(round(linec(stoppos,1)),round(linec(stoppos,2))) = 1;
%        end
        F = bwdist(mask) >= 0.5*circle.*diamo;
        F = imcomplement(F);
    end
    
    if (split == 1)
        if (circle > 0)
            stoppos = length(linec); stopc1 = length(total1); stopc2 = length(total2);
        end
        roi1 = vertcat(total1(startc1:stopc1,:), linec(stoppos:-1:startpos,:));
        roi2 = vertcat(total2(startc2:stopc2,:), linec(stoppos:-1:startpos,:));
%        if (starti < distc_t*pixelsize)
%            roi1 = vertcat(boundb(range2(1):postotal1(2),:),roi1,boundb(range2(1),:));
%            roi2 = vertcat(boundb(range2(1):-1:postotal2(1),:),roi2,boundb(range2(1),:));
%        end
        FS1 = F.*poly2mask(roi1(:,2),roi1(:,1),Esize(1),Esize(2));
        FS2 = F.*poly2mask(roi2(:,2),roi2(:,1),Esize(1),Esize(2));
    end
end