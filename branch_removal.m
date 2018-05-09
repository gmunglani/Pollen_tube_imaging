function [Q2, Qe,Qarea] = branch_removal(Q,Qbf,Qel,angdiff,close_dist)

Q2 = Q;
Qbr = Qbf(:,1);
Qbc = Qbf(:,2);
Qarea = 1;

while(length(Qbr) > 0)
    Qtmp = Q2;
    Q2(Qbr(1)-1:Qbr(1)+1,Qbc(1)-1:Qbc(1)+1) = 0;
    Q2label = bwlabel(Q2);
    Q2stat = regionprops(Q2label,'Area','PixelList','Centroid');   
    Q2area = []; Q2cen = [];
    for i = 1:length(Q2stat)
        Q2area = [Q2area Q2stat(i).Area];
        Q2cen = [Q2cen Q2stat(i).Centroid(1)];
    end

    if (angdiff > 0)
        [tmp, Q2ai] = sort(Q2area);
        Qarea = bwarea(Q2label(Q2label == Q2ai(1)))/bwarea(Q2label(Q2label == Q2ai(2)));
    end
    
    Qea = [];
    if (length(Q2stat) < 3)
        [tmp Q2pos] = min(Q2cen);
        Qea = intersect(circshift(Qel,1,2),Q2stat(Q2pos).PixelList,'rows');
        
        Q2(Qbr(1),Qbc(1)) = 1;
        Q2 = bwmorph(Q2,'bridge');
        Q2 = bwmorph(Q2,'clean');
        Q2 = bwmorph(Q2,'thin');
        Qbr(1) = []; Qbc(1) = [];
    else
        [tmp, Q2ci] = sort(Q2cen);
        Q2a1 = Q2stat(Q2ci(1)).PixelList;
        Q2a2 = Q2stat(Q2ci(2)).PixelList;
        Qea = [intersect(circshift(Qel,1,2),Q2a1,'rows'); intersect(circshift(Qel,1,2),Q2a2,'rows')];
        
        Q2m = zeros(size(Q2)); Q2m = createcircles(Q2m,[Qbc(1) Qbr(1)],20);
        Q3label = Q2m.*Q2label;
        Q3stat = regionprops(Q3label,'Orientation','Area');
        
        Q3area = [];
        for i = 1:length(Q3stat)
            Q3area = [Q3area Q3stat(i).Area];
        end
        
        [Q3areamin Q3areai] = min(Q3area);
        
        if (Q3areamin < 4) 
            Q3min = Q3areai;
            Q3diffs = [1 2 1];
        else
            % Angle check
            Q3ang = [Q3stat(1).Orientation Q3stat(2).Orientation Q3stat(3).Orientation];
            Q3diff = diff(Q3ang); Q3diff(3) = Q3ang(1)-Q3ang(3);
            u = Q3diff >= 90;
            Q3diff(u) = 180 - Q3diff(u);
            l = Q3diff <= -90;
            Q3diff(l) = 180 + Q3diff(l);
            Q3diff = abs(circshift(Q3diff,2));
            [Q3diffs, Q3di] = sort(Q3diff);
            if (close_dist == 2 && Q2ci(1) == Q3di(1))
                Q3cdiff = diff(Q2cen); Q3cdiff(3) = Q2cen(1)-Q2cen(3);
                [tmp,Q3cdiffp] = min(abs(Q3cdiff));
                
                Q3min = Q3di(2);
            else
                Q3min =  Q3di(1);
            end
        end
        
        % Find edges
        Q2ci(find(Q2ci == Q3min)) = [];
        Q3a1 = Q2stat(Q2ci(1)).PixelList;
        Q3a2 = Q2stat(Q2ci(2)).PixelList;
        
        % Compare areas of 2 main branches
        Qeb1 = intersect(Qea,Q3a1,'rows');
        Qeb2 = intersect(Qea,Q3a2,'rows');
        [tmp,Qep1] = intersect(Qea,Qeb1,'rows');
        [tmp,Qep2] = intersect(Qea,Qeb2,'rows');
        if (Qep1 == 2)
            Qea = circshift(Qea,1,1);
            Qarea = 1/Qarea;
        elseif (Qep2 == 2)
            Qea = circshift(Qea,1,1);
            Qarea = 1/Qarea;
        end
        
        if(Q3diffs(2) - Q3diffs(1) >= angdiff*0.25 && Q3diffs(3) > angdiff && close_dist > 0)
            Q3label(find(Q3label == Q3min)) = 0;
            Q3label(find(Q3label > 0 )) = 1;
            Q4bridge = Q3label;
            Q4bridge(Qbr(1),Qbc(1)) = 1;
            Q4bridge = bwmorph(Q4bridge,'bridge');
            Q4bridge = bwmorph(Q4bridge,'clean');
            Q4bridge = bwmorph(Q4bridge,'thin');
            Q4extra = Q4bridge - Q3label;
            Q2label(find(Q2label == Q3min)) = 0;
            Q2 = Q2label + Q4extra;
            Q2(find(Q2 > 0 )) = 1;
            if (size(Qea,1) > 1) Qea(2,:) = []; end
        else
            Q2 = Qtmp;         
        end
        Qbr(1) = []; Qbc(1) = [];
    end
end
if (angdiff == 0)
    Qei = bwmorph(Q2,'endpoints');
    [Qer,Qec] = find(Qei > 0);
    Qe = [Qer Qec];
    [tmp,Qepos] = max(Qec);
    Qe(Qepos,:) = [];
else    
    Qe = circshift(Qea,1,2);
end