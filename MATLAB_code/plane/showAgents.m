function output = showAgents(X,U,tar,DT,fixframe,Ndogs,frameinrad,t,equilibrium,alphaHull,timeStratifiedX,plotTSX,expDecay,colours,pallette,scalarField,subflocksCell)


purple = [128, 0 , 128] / 255;

clipToVoronoi=false;

% NEED CM DEFINED HERE!

fig = figure(1);
scatter(X(1:Ndogs,1), X(1:Ndogs,2), 'filled','red')
hold on
N = numel(X(:,1));


%---------------------------%
%for rainbow sheep uncomment
%scatter(X(Ndogs+1:N,1),X(Ndogs+1:N,2), 36,colours,'filled')
%scatter(X(Ndogs+1:N,1),X(Ndogs+1:N,2),10,"blue",'filled')
% disp(numel(subflocksCell))
for i=1:numel(subflocksCell)
    scatter(subflocksCell{i}(:,1),subflocksCell{i}(:,2),'filled')
end

%---------------------------%
if plotTSX
    timeStratifiedX( find(~any(timeStratifiedX')), : ) = [];
    mem = scatter(timeStratifiedX(:,1),timeStratifiedX(:,2),40, purple,'filled',"s");
    % mem.AlphaData = 5.*expDecay;
    % mem.MarkerFaceAlpha = 'flat';
    % CM = sum(expDecay.*timeStratifiedX,1)./sum(expDecay);
    % scatter(CM(1),CM(2),'red','s')

end
%---------------------------%

colormap(pallette);

hold on
% scatter(positions(:,1), positions(:,2), 'green')
% scatter(equil(:,1), equil(:,2))
% scatter(CM(:,1), CM(:,2),'black')

if Ndogs>0
    scatter(equilibrium(1), equilibrium(2), 'green')
end
% shp = alphaShape(Xsheep,sqrt(N));
% [~,P] = boundaryFacets(shp);
% plot(shp)
%plot(P(:,1),P(:,2))
% plot(alphaHull(:,1),alphaHull(:,2),'Color',purple)

%---------------------------%
qs = 3;

quiver(X(:,1),X(:,2),qs*U(:,1),qs*U(:,2),'off', 'k');

C = tar.Components;
for m=1:size(C,1)
    if isa(C{m},'polyshape')
        plot(C{m},'EdgeColor','none');
    end
end



clip = [min(X,[],1); max(X,[],1)];
if clipToVoronoi==true
    V = voronoiDiagram(DT);
    clip = [min([V; X],[],1); max([V(2:end,:); X],[],1)];
end
clip_width = 1.2*max(clip(2,:)-clip(1,:));
clip_center = mean(clip,1);
clip_lim(1:4) = clip_center(ceil((1:4)/2)) - rescale(mod(1:4,2),-1,1)*.5*clip_width;

axis equal

frame = [-frameinrad frameinrad -frameinrad frameinrad];
if ~fixframe
    com = mean(X);
    rmed = sqrt(median((X(:,1)-com(1)).^2+(X(:,2)-com(2)).^2));
    frame = [com(1)-3*rmed com(1)+3*rmed com(2)-3*rmed com(2)+3*rmed];
end
axis(frame);

%---------------------------%
%preference field
lim = axis;
hold on
if ~strcmp(scalarField,'zero')
    lim = axis;
    ax = linspace(lim(1),lim(2));
    [X,Y,Z] = getMeshScalarField(ax,scalarField);
    levels = linspace(min(Z,[],"all"),max(Z,[],"all"),20);
    levels(end+1) = 0;
    if (strcmp(scalarField,'fence'))
        yFence = linspace(lim(1),15,100);
        scatter(yFence,20.*ones(size(yFence)),20,'filled','s','black')
        xFence = linspace(lim(3),15,100);
        scatter(20.*ones(size(xFence)),xFence,20,'filled','s','black')

        %contour(X,Y,Z,6,'black','LineWidth',0.2)
    elseif (strcmp(scalarField,'fenceNoGap'))
        yFence = linspace(lim(1),20,100);
        scatter(yFence,20.*ones(size(yFence)),20,'filled','s','black')
        xFence = linspace(lim(3),20,100);
        scatter(20.*ones(size(xFence)),xFence,20,'filled','s','black')
    elseif (strcmp(scalarField,'infiniteFence'))
        xFence = linspace(lim(3),lim(4),100);
        scatter(20.*ones(size(xFence)),xFence,20,'filled','s','black')
    else
        contour(X,Y,Z,levels,'black','LineWidth',1.5)

    end
    axis(lim)
end
%---------------------------%

title(['$t=',num2str(t),'$'], 'Interpreter', 'Latex');
hold off

output = 0;