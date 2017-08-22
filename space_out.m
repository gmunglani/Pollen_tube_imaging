function xx1 = space_out(xx1)

rept1 = diff([0 find(diff(xx1)) numel(xx1)]);
repos1 = find(rept1>1);
for d = 1:length(repos1)
    if (repos1(d) == 1)
        repfi1 = linspace(xx1(repos1(d))-1,xx1(repos1(d)+rept1(repos1(d))),rept1(repos1(d))+2);
        xx1(repos1(d):repos1(d)+rept1(repos1(d))) = repfi1(2:end);
    else
        repfi1 = linspace(xx1(repos1(d)-1),xx1(repos1(d)+rept1(repos1(d))),rept1(repos1(d))+2);
        xx1(repos1(d)-1:repos1(d)+rept1(repos1(d))) = repfi1;
    end
end