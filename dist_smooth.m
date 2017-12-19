function [avgxy,diffxy] = dist_smooth(ref,loc,bound,avg,avg_diff)

leng = []; 
j = 1;
for i = loc(1):loc(end)   
    leng(j) = pdist2(bound(i,:),ref);
    j=j+1;
end
avgxy = movmean(leng,avg);
diffxy = movmean(diff(avgxy),avg_diff);

