% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
%% Code for reproducing results in Figure 3
p = 50; % Change to p=2,10,1700 for other results in figure
% Load and modify data
filename='bild2.png';
[graph,values,m,n,image]=read_Image(filename);
g=image(:);
U_sq=zeros(m,n);
U_sq((round(m/2)-19):(round(m/2)+20),(round(n/2)-19):(round(n/2)+20))=1;
U_sq=logical(U_sq);
U=~reshape(U_sq,m*n,1)';
image2=image;
image2(repmat(U_sq,1,1,3))=0;
imagesc(image);
title('Initial image')
axis equal
axis off

%% Compute componentwise tight extension
mask = ~U';
erg=reshape(g,3,size(g,1)/3)';
erg(:,1) = reshape(image(:,:,1),1,[]);
erg(:,2) = reshape(image(:,:,2),1,[]);
erg(:,3) = reshape(image(:,:,3),1,[]);
erg2=ComponentwiseTight(graph,erg,mask,'epsilon',1e-9,...
    'maxIterations',10000);
Llex_com = computeL_lex(graph,~mask,erg2);
im_out=reshape(erg2',[],1);
write_Image('Square_RGB_Comp.png',im_out,m,n,1);
title('Componentwise tight extension')

%% Compute iterated midrange filter
img = erg';
img = img(:);
for i =1:500
    img = refine_solution(graph,U,img,3);
end
im_out=reshape(img',[],1);
write_Image('Square_RGB_Midr.png',im_out,m,n,1);
title('Iterated midrange filter')

%% Compute approximation of tight extension
[erg,Llex]=ApproximateTight(graph,erg,mask,p);
im_out=reshape(erg',[],1);
write_Image('Square_RGB_Tight.png',im_out,m,n,1);
title(['Approximation of tight extension with p=',num2str(p)])
