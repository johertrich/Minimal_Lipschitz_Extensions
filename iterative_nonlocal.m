% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
% Algorithm for image inpainting
% INPUT:
% img       - image for inpainting
% mask      - mask indicating missing pixels (1 corresponds to unknwown)
% OPTIONAL:
% epsilon           - tolerance
% maxIterations     - maximal iteration number
% laplacian         - functionhandle for Laplaceoperator (see correponding
%                     functions)
% patch_p           - radius (each direction) of patches which are compared
% patch_r           - radius (each direction) in which patches are compared
% nNeighbors        - desired number of neighbours
% sigma             - weighting parameter, smaller corresponds to less
%                     focus on deviations in the weights
% local_weight      - weighting to enforce some localness of neighbours
% plot              - 1 if results should be plotted after each step
% min_Steps         - minimum number of iterations
% p_max             - p_max is the used approximation
% mode              - 0 default mode
%                     1 if pixels should remain fixed after computing them
% OUTPUT:
% out               - inpainted image
function out=iterative_nonlocal(img,mask,varargin)
[size_m,size_n,dimensions]=size(img);
p=inputParser;
addOptional(p,'epsilon',1e-5);
addOptional(p,'maxIterations',10000);
addOptional(p,'laplacian',@(g,f)infLaplacian(g,f,0.5));
addOptional(p,'patch_p',7);
addOptional(p,'patch_r',15);
addOptional(p,'nNeighbors',8);
addOptional(p,'sigma',0.2);
addOptional(p,'local_weight',0);
addOptional(p,'plot',0);
addOptional(p,'graph',graphGen4_neighbours(size_m,size_n));
addOptional(p,'min_Steps',30);
addOptional(p,'p_max',200);
addOptional(p,'mode',0);
parse(p,varargin{:});
s=p.Results;

if s.plot==1
    figure('units','normalized','outerposition',[0 0 1 1])
end

% Reshape everything to vector form
f=zeros(size_m*size_n,dimensions);
for i=1:dimensions
    f(:,i)=reshape(img(:,:,i),[],1);
end
step=1;
if s.mode==1
    mask2=mask(:);
    break_condition = max(mask2)>0;
else
    mask2=zeros(size_m*size_n,1);
    break_condition = (max(mask2)>0 || step-1<s.min_Steps) && (step <10);
end

% Iterate until all pixels are known
while break_condition
    disp(['step ' num2str(step) ' left ' num2str(sum(mask2)) ' elements']);
    
    % compute nonlocal graph
    aeusseres=(s.graph*(1-mask2))~=0; % pixels with known neighbour in graph
    rand_maske=aeusseres.*mask2; % unknown pixels with neighbour in graph
    if s.mode==1
        graph2=graphGen_nonlocal(img,reshape(mask2,size_m,size_n),...
            'patchSize',s.patch_p,'radius',s.patch_r,'nNeighbors',s.nNeighbors,...
            'sigma',s.sigma,'needed_area',logical(reshape(rand_maske,size_m,size_n)));
    else
        graph2=graphGen_nonlocal(img,reshape(mask2,size_m,size_n),...
            'patchSize',s.patch_p,'radius',s.patch_r,'nNeighbors',s.nNeighbors,...
            'sigma',s.sigma,'needed_area',mask);
    end
    
    % add relevant part of (local) graph to non-local one
    anz_elem_ausserhalb_rand_maske=sum(1-rand_maske);
    graph_auf_rand_maske_adj=s.graph;
    graph_auf_rand_maske_adj(logical(1-rand_maske),logical(1-rand_maske))=...
        sparse(anz_elem_ausserhalb_rand_maske,anz_elem_ausserhalb_rand_maske);
    graph2=graph2 + s.local_weight*graph_auf_rand_maske_adj;
    
    disp('graph generated');
    % Compute tight extension with p_max and update break_condition
    if s.mode==1
        f=ApproximateTight(graph2,f,mask2,s.p_max);
        mask2=(mask2-((graph2*(1-mask2))~=0).*((s.graph*(1-mask2)~=0)))>0;
        break_condition = max(mask2)>0;
    else
        f=ApproximateTight(graph2,f,mask(:),s.p_max);
        mask2=(mask2-((graph2*(1-mask2))~=0))>0;
        break_condition = (max(mask2)>0 || step-1<s.min_Steps) && (step <10);
    end

    img=zeros(size_m,size_n,dimensions);
    for i=1:dimensions
        img(:,:,i)=reshape(f(:,i),[size_m,size_n]);
    end
    % show result
    if s.plot==1
        imagesc(img);
        axis image
        axis off
        title(['After step ' num2str(step)])
    end
    drawnow
    step=step+1;
end
disp(['Terminated after ' num2str(step) ' steps']);
out=f;
end

% Generates local neighbourhood graph
function out=graphGen4_neighbours(m,n)
    Bm=spdiags([ones(m,1) , zeros(m,1) , ones(m,1)], [-1 0 1], m, m);
    Bn=spdiags([ones(n,1) , zeros(n,1) , ones(n,1)], [-1 0 1], n, n);
    out=kron(speye(n),Bm)+kron(Bn,speye(m));
end
