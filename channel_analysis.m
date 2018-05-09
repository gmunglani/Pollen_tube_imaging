function [Mmin Mmax Mmin_prc Mmax_prc Mmed] = channel_analysis(M,smp)

for count = 1:smp
    C = M(:,:,count);
    Mmax(count) = max(max(C));
    Carray = reshape(C,size(C,2)*size(C,1),1);
    Carray = Carray(Carray>0);
    Mmax_prc(count) = prctile(Carray,95);
    Mmed(count) = median(Carray);
    if (sum(C(C>0)) > 1) Mmin(count) = min(C(C>0));
    else Mmin(count) = 0;
    end
    Mmin_prc(count) = prctile(Carray,5);
end

