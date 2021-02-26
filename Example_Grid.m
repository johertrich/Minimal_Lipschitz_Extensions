% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
%% Code for reproducing images in Figure 2.
p = 200; % Change to p=2,50,2400 for other results in figure
%Load Data
load('DataGitter.mat')
mask=~U';

%% Compute componentwise tight extension
erg=reshape(g,2,size(g,1)/2)';
erg_component=ComponentwiseTight(graph,erg,mask,'epsilon',1e-10,...
    'maxIterations',10000);
% Plot result
g_out=erg_component';
plot(g_out(1,U),g_out(2,U),'o','MarkerSize',10);
hold on;
plot(g_out(1,~U),g_out(2,~U),'x','MarkerSize',14);
for i=1:size(graph,1)
    for j=find(graph(i,:))
        plot(g_out(1,[i,j]),g_out(2,[i,j]),'b-','LineWidth',1.2);
    end
end
title('Componentwise tight extension')
axis equal;
hold off;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%% Compute iterated midrange filter
img = erg';
img = img(:);
for i =1:500
    img = refine_solution(graph,U,img,2);
end
erg_midrange = reshape(img,2,[])';
% Plot result
figure;
g_out=erg_midrange';
plot(g_out(1,U),g_out(2,U),'o','MarkerSize',10);
hold on;
plot(g_out(1,~U),g_out(2,~U),'x','MarkerSize',14);
for i=1:size(graph,1)
    for j=find(graph(i,:))
        plot(g_out(1,[i,j]),g_out(2,[i,j]),'b-','LineWidth',1.2);
    end
end
title('Iterated midrange filter')
axis equal;
hold off;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% Compute approximation of tight extension
erg_tight=ApproximateTight(graph,erg,mask,p);
% Plot result
figure;
g_out = erg_tight';
plot(g_out(1,U),g_out(2,U),'o','MarkerSize',10);
hold on;
plot(g_out(1,~U),g_out(2,~U),'x','MarkerSize',14);
for i=1:size(graph,1)
    for j=find(graph(i,:))
        plot(g_out(1,[i,j]),g_out(2,[i,j]),'b-','LineWidth',1.2);
    end
end
title(['Approximation of tight extension with p=',num2str(p)])
axis equal;
hold off;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
