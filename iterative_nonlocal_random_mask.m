% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
% Algorithm for image inpainting with random masks
% INPUT:
% img       - image for inpainting
% mask      - mask indicating missing pixels (1 corresponds to unknwown)
% OPTIONAL:
% delta_t           - stepsize for diffusion equation
% epsilon           - tolerance
% maxIterations     - maximal iteration number
% laplacian         - functionhandle for Laplaceoperator (see correponding
%                     functions)
% patch_p           - radius (each direction) of patches which are compared
% patch_r           - radius (each direction) in which patches are compared
% nNeighbors        - desired number of neighbours
% local_weight      - weighting to enforce some localness of neighbours
% plot              - 1 if results should be plotted after each step
% min_Steps         - minimum number of iterations
% p_max             - p_max is the used approximation
% OUTPUT:
% out               - inpainted image
function out=iterative_nonlocal_random_mask(img,mask,size_m,size_n,varargin)
[~,dimensions]=size(img);
p=inputParser;
addOptional(p,'delta_t',0.6);
addOptional(p,'epsilon',1e-5);
addOptional(p,'maxIterations',10000);
addOptional(p,'laplacian',@(g,f)infLaplacian(g,f,0.5));
addOptional(p,'patch_p',11);
addOptional(p,'patch_r',15);
addOptional(p,'nNeighbors',8);
addOptional(p,'local_weight',0.5);
addOptional(p,'plot',0);
addOptional(p,'min_Steps',15);
addOptional(p,'p_max',200);
parse(p,varargin{:});
s=p.Results;
s.graph=graphGen4_neighbours(size_m,size_n);

if s.plot==1
    figure('units','normalized','outerposition',[0 0 1 1])
end
step=1;
mask2 = mask;
% Iterate until all pixels are known
while (max(mask)>0 || step-1<s.min_Steps) && (step-1<15)
    % compute nonlocal graph
    graph2=graphGen_mex(img,mask2,size_m,size_n,'patchSize',s.patch_p,...
        'radius',s.patch_r,'nNeighbors',s.nNeighbors);
    aeusseres=(s.graph*(1-mask2))~=0; % pixels with known neighbour in graph
    rand_maske=aeusseres.*mask2; % unknown pixels with neighbour in graph
    
    % add relevant part of (local) graph to non-local one
    anz_elem_ausserhalb_rand_maske=sum(1-rand_maske);
    graph_auf_rand_maske_adj=s.graph;
    graph_auf_rand_maske_adj(logical(1-rand_maske),logical(1-rand_maske))=...
        sparse(anz_elem_ausserhalb_rand_maske,anz_elem_ausserhalb_rand_maske);
    graph2=graph2+s.local_weight*graph_auf_rand_maske_adj;
    s.local_weight = max(s.local_weight*0.8,0.05);
    disp('Graph generated');
    
    % Other extension could be inserted here
    img=ApproximateTight(sqrt(graph2),img,mask2,s.p_max);
    
    % Plot result
    if s.plot==1
        img2=zeros(size_m,size_n,dimensions);
        for i=1:dimensions
            img2(:,:,i)=reshape(img(:,i),[size_m,size_n]);
        end
        imshow(img2,[]);
        axis image
        axis off
        title(['After step ' num2str(step)])
    end
    drawnow
    step=step+1;
end
disp(['Terminated after ' num2str(step) ' steps']);
out=img;
end

% Generates local neighbourhood graph
function out=graphGen4_neighbours(m,n)
    Bm=spdiags([ones(m,1) , zeros(m,1) , ones(m,1)], [-1 0 1], m, m);
    Bn=spdiags([ones(n,1) , zeros(n,1) , ones(n,1)], [-1 0 1], n, n);
    out=kron(speye(n),Bm)+kron(Bn,speye(m));
end
