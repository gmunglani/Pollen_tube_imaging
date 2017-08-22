function [totals1 totals2] = angle_correction(totals1,totals2,maxy,diam)

reg1 = find(totals1(:,2)>(maxy-diam));
reg2 = find(totals2(:,2)>(maxy-diam));

[cpfit1,err] = polyfit(totals1(reg1(1):reg1(end),2),totals1(reg1(1):reg1(end),1),1);
[cpfit2,err] = polyfit(totals2(reg2(1):reg2(end),2),totals2(reg2(1):reg2(end),1),1);
    
cpfit = (cpfit1(1)+cpfit2(1))*0.5;
deg = pi/2 - atan(abs(cpfit));
xsub = cos(deg(1))*diam*sin(deg);
    
if (cpfit > 0) totals2(find(totals2(:,2) > (maxy-xsub)),:) = [];
elseif (cpfit < 0) totals1(find(totals1(:,2) > (maxy-xsub)),:) = [];
end