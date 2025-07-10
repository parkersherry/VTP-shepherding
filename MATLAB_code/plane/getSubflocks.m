function S = getSubflocks(X, XAtCurrSheepnbhd, dogL, N, Ndogs, subflocksCell, DT_on)

% If we use vertices, https://www.mathworks.com/help/matlab/ref/triangulation.nearestneighbor.html

arguments
    X
    XAtCurrSheepnbhd
    dogL
    N
    Ndogs
    subflocksCell
    DT_on
end
subflockA = {};
subflockB = {};
% dbstack

thresh = 9;
% [A,B] = deal([]);
save = 0;

longestLength = 0;
LL_index = 0;

DT = delaunayTriangulation(XAtCurrSheepnbhd);
E = edges(DT);
numEdges = (numel(E)/2);
for i = 1:numEdges
    startpoint = XAtCurrSheepnbhd(E(i,1),:);
    endpoint = XAtCurrSheepnbhd(E(i,2),:);
    length = vecnorm(startpoint - endpoint);
    if length > longestLength
        longestLength = length;
        LL_index = i;
        a = startpoint;
        b = endpoint;
    end
end
if longestLength > thresh
    % [a,b] = deal(X(E(LL_index,1),:), X(E(LL_index,2),:));
    sort_herd = sign(vecnorm(X-a,2,2) - vecnorm(X-b,2,2));
    [A,B] = deal(find(sort_herd < 0), find(sort_herd > 0));

    %find shortest length between CH_A and CH_B, called Gamma
    [X_A, X_B] = deal(X(A,:), X(B,:));
    % fig = figure(1);
    % hold on
    % scatter(X_A(:,1),X_A(:,2),'filled','red')
    % scatter(X_B(:,1),X_B(:,2),'filled','blue')

    [CH_A, CH_B] = deal(convhull(X_A), convhull(X_B));
    [tarA, tarB] = deal(Target(X_A(CH_A,:)), Target(X_B(CH_B,:)));
    [minIndexA, minIndexB] = deal(0);
    minLength = Inf;

    for j = 1:numel(CH_A)
        h = homeToTarget(tarB, X_A(CH_A(j),:));
        if vecnorm(h) < minLength
            minLength = vecnorm(h);
            minIndexA = j;
        end
    end
    minLength = Inf;
    for j = 1:numel(CH_B)
        h = homeToTarget(tarA, X_B(CH_B(j),:));
        if vecnorm(h) < minLength
            minLength = vecnorm(h);
            minIndexB = j;
        end
    end

    nbrIndexMinusA = minIndexA-1;
    nbrIndexPlusA = minIndexA+1;
    nbrIndexMinusB = minIndexB-1;
    nbrIndexPlusB = minIndexB+1;

    if nbrIndexMinusA == 0
        nbrIndexMinusA = numel(CH_A);
    end

    if nbrIndexPlusA > numel(CH_A)
        nbrIndexPlusA = 1;
    end

    if nbrIndexMinusB == 0
        nbrIndexMinusB = numel(CH_B);
    end

    if nbrIndexPlusB > numel(CH_B)
        nbrIndexPlusB = 1;
    end

    [P_A, P_B, nA1, nA2, nB1, nB2] = deal(CH_A(minIndexA), CH_B(minIndexB), CH_A(nbrIndexMinusA), CH_A(nbrIndexPlusA), CH_B(nbrIndexMinusB), CH_B(nbrIndexPlusB));

    Gamma1 = lineDist(X_A, X_B, P_A, P_B, nA1, nB1);
    Gamma2 = lineDist(X_A, X_B, P_A, P_B, nA1, nB2);
    Gamma3 = lineDist(X_A, X_B, P_A, P_B, nA2, nB1);
    Gamma4 = lineDist(X_A, X_B, P_A, P_B, nA2, nB2);

    Gamma = min([Gamma1 Gamma2 Gamma3 Gamma4]);


    % if DT_on == 1
    %     fig1 = figure(1);
    %     scatter(X_A(:,1), X_A(:,2), 'filled')
    %     hold on
    %     scatter(X_B(:,1), X_B(:,2), 'filled')
    %     % triplot(DT)
    %     scatter(a(1), a(2), 's')
    %     scatter(b(1), b(2), 's')
    % end

    disp(Gamma)

    if Gamma > dogL
        save = 1;
        % do this again until there are no more subflocks, i.e.
        ACurr = intersect(X_A,XAtCurrSheepnbhd, 'rows');
        subflockA = getSubflocks(X_A,ACurr,dogL,N,Ndogs, subflocksCell,1);
        BCurr = intersect(X_B,XAtCurrSheepnbhd, 'rows');
        subflockB = getSubflocks(X_B,BCurr,dogL,N,Ndogs,subflocksCell,1);
    end
end
if (save==0)
    subflocksCell{end+1} = X;
end

S = [subflocksCell;subflockA;subflockB];