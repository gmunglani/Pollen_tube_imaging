function [poscross1, poscross2, distcf] = line_continuity(poscross1, poscross2, type, distc)

if (type == 1)
    poscrossa = poscross1;
    poscrossb = poscross2;
else
    poscrossb = poscross1;
    poscrossa = poscross2;
end

if (type == 1) distc = horzcat(0,distc); end

dposcrossa = find(diff(poscrossa) <= 0);
while ~isempty(dposcrossa) 
    poscrossa(dposcrossa+1) = [];
    poscrossb(dposcrossa+1) = [];
    distc(dposcrossa+1) = []; 
    dposcrossa = find(diff(poscrossa) <= 0);
end

% dposcrossao = find(diff(poscrossa) < 0)
% while ~isempty(dposcrossao)
%     for a = 1:length(dposcrossao)
%         if (dposcrossao(a) > 1)
%             posnew1 = round((poscrossa(dposcrossao(a)+1) + poscrossa(dposcrossao(a)-1))*0.5);
%             if (posnew1 == poscrossa(dposcrossao(a)+1) | posnew1 == poscrossa(dposcrossao(a)-1))
%                 tmp1 = poscrossa(dposcrossao(a));
%                 poscrossa(dposcrossao(a)) = poscrossa(dposcrossao(a)+1);
%                 poscrossa(dposcrossao(a)+1) = tmp1;
%             else
%                 poscrossa(dposcrossao(a)) = posnew1;
%             end
%         else
%             poscrossa(dposcrossao(a)) = poscrossa(dposcrossao(a)+1);
%         end
%     end
%end

if (type == 1)
    poscross1 = poscrossa;
    poscross2 = poscrossb;
else
    poscross2 = poscrossa;
    poscross1 = poscrossb;
end
distcf = distc;