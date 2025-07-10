function D = lineDist(XA, XB, PA, PB, nA, nB)

arguments
    XA
    XB
    PA
    PB
    nA
    nB
end

[clPtA, clPtB] = deal(XA(PA,:), XB(PB,:));
[linePtA, linePtB] = deal(XA(nA,:), XB(nB,:));

if clPtA(1) == linePtA(1) && clPtB(1) ~= linePtB(1)
    mB =  (clPtB(2) - linePtB(2))/(clPtB(1) - linePtB(1));
    interceptB = clPtB(2)-mB*clPtB(1);
    lineIntersect = [clPtA(1) mB*clPtA(1)+interceptB];

    endPtsA = [clPtA; linePtA];
    endPtsB = [clPtB; linePtB];

    [~,endPA]= min([vecnorm(lineIntersect - clPtA) vecnorm(lineIntersect - linePtA)]);
    [~,endPB]= min([vecnorm(lineIntersect - clPtB) vecnorm(lineIntersect - linePtB)]);

    if (endPtsA(1,2) > endPtsB(endPB,2) && endPtsB(endPB,2) > endPtsA(2,2)) || (endPtsA(1,2) < endPtsB(endPB,2) && endPtsB(endPB,2) < endPtsA(2,2))
        D = abs(endPtsB(endPB,1) - endPtsA(1,1));
    else
        D = abs(mB*endPtsA(endPA,1)-endPtsA(endPA,2)+interceptB)/sqrt(mB^2+1);
    end

elseif clPtA(1) ~= linePtA(1) && clPtB(1) == linePtB(1)
    mA =  (clPtA(2) - linePtA(2))/(clPtA(1) - linePtA(1));
    interceptA = clPtA(2)-mA*clPtA(1);
    lineIntersect = [clPtB(1) mA*clPtB(1)+interceptA];

    endPtsA = [clPtA; linePtA];
    endPtsB = [clPtB; linePtB];

    [~,endPA]= min([vecnorm(lineIntersect - clPtA) vecnorm(lineIntersect - linePtA)]);
    [~,endPB]= min([vecnorm(lineIntersect - clPtB) vecnorm(lineIntersect - linePtB)]);

    if (endPtsB(1,2) > endPtsA(endPA,2) && endPtsA(endPA,2) > endPtsB(2,2)) || (endPtsB(1,2) < endPtsA(endPA,2) && endPtsA(endPA,2) < endPtsB(2,2))
        D = abs(endPtsA(endPA,1) - endPtsB(1,1));
    else
        D = abs(mA*endPtsB(endPB,1)-endPtsB(endPB,2)+interceptA)/sqrt(mA^2+1);
    end

elseif clPtA(1) == linePtA(1) && clPtB(1) == linePtB(1)
    lineArrA = sort([clPtA; linePtA], 1);
    lineArrB = sort([clPtB; linePtB], 1);

    y1 = lineArrA(1,2);
    y2 = lineArrA(2,2);
    y3 = lineArrB(1,2);
    y4 = lineArrB(2,2);

    if (y1 <= y3 && y3 <= y2) || (y1 <= y4 && y4 <= y2) || (y3<=y1 && y2<=y4) || (y1<=y3 && y4<=y2)
        D = abs(clPtA(1) - clPtB(1));
    elseif y1 > y4
        D = vecnorm(lineArrA(1,:) - lineArrB(2,:));
    else
        D = vecnorm(lineArrA(2,:) - lineArrB(1,:));
    end


else
    mA =  (clPtA(2) - linePtA(2))/(clPtA(1) - linePtA(1));
    mB =  (clPtB(2) - linePtB(2))/(clPtB(1) - linePtB(1));
    interceptA = clPtA(2)-mA*clPtA(1);
    interceptB = clPtB(2)-mB*clPtB(1);

    if mA == mB
        tempVec = clPtA - linePtA;
        theta = atan(tempVec(1)/tempVec(2));
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

        clPtA = (R*clPtA')';
        linePtA = (R*linePtA')';
        clPtB = (R*clPtB')';
        linePtB = (R*linePtB')';

        lineArrA = sort([clPtA; linePtA], 1);
        lineArrB = sort([clPtB; linePtB], 1);

        y1 = lineArrA(1,2);
        y2 = lineArrA(2,2);
        y3 = lineArrB(1,2);
        y4 = lineArrB(2,2);

        if (y1 <= y3 && y3 <= y2) || (y1 <= y4 && y4 <= y2) || (y3<=y1 && y2<=y4) || (y1<=y3 && y4<=y2)
            D = abs(clPtA(1) - clPtB(1));
        elseif y1 > y4
            D = vecnorm(lineArrA(1,:) - lineArrB(2,:));
        else
            D = vecnorm(lineArrA(2,:) - lineArrB(1,:));
        end
    else

        lineIntersect = [(interceptA - interceptB)/(mB - mA) mB*((interceptA - interceptB)/(mB - mA))+interceptB];

        endPtsA = [clPtA; linePtA];
        endPtsB = [clPtB; linePtB];

        [~,endPA]= min([vecnorm(lineIntersect - clPtA) vecnorm(lineIntersect - linePtA)]);
        [~,endPB]= min([vecnorm(lineIntersect - clPtB) vecnorm(lineIntersect - linePtB)]);

        D1 = abs(mB*endPtsA(endPA,1)-endPtsA(endPA,2)+interceptB)/sqrt(mB^2+1);
        D2 = abs(mA*endPtsB(endPB,1)-endPtsB(endPB,2)+interceptA)/sqrt(mA^2+1);

        if D1==0
            D=D2;
        elseif D2==0
            D=D1;
        else
            D = min([D1 D2]);
        end
    end
end

% get endpoint that is closer to lineIntersect


