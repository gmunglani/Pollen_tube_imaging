keep00 = [];keep11 = [];
keep10 = [];keep01 = [];

pos00 = [];pos11 = [];
pos10 = [];pos01 = [];

for a=1:count
    if (fix(a) == 2 && cacl(a) > 1)
        keep11 = [keep11; phin(a)];
        pos11 = [pos11; a];
    elseif (fix(a) == 1 && cacl(a) >1)
        keep10 = [keep10; phin(a)];
        pos10 = [pos10; a];
    elseif (fix(a) == 2 && cacl(a) <1)
        keep01 = [keep01; phin(a)];
        pos01 = [pos01; a];
    else
        keep00 = [keep00; phin(a)];
        pos00 = [pos00; a];
    end
end

figure
plot(pos00,keep00,'b*')
hold on
plot(pos10,keep10,'k*')
plot(pos01,keep01,'g*')
plot(pos11,keep11,'r*')