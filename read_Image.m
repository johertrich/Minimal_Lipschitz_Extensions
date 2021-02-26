% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
% Loads an RGB-image, creates an 4-neighbourhood graph on it and
% converts the values of the image into a vectorized form.
% INPUT:
%   filename - path to the image
% OUTPUT:
%   graph - adjacency matrix of a 4-neighbourhood graph
%   values - contains the RGB image values in a vectorized form
%   m,n - size of the image
function [graph,values,m,n,image]=read_Image(filename)
    im=imread(filename);
    r=double(im(:,:,1))/256;
    g=double(im(:,:,2))/256;
    b=double(im(:,:,3))/256;
    image(:,:,1)=r;
    image(:,:,2)=g;
    image(:,:,3)=b;
    
    [m,n]=size(r);
    
    rgb_vector=[reshape(r,m*n,1) reshape(g,m*n,1), reshape(b,m*n,1)]';
    
    values=reshape(rgb_vector,3*m*n,1);
    
    graph=graphGen4_neighbours(size(r,1),size(r,2));
end

% Generates local neighbourhood graph
function out=graphGen4_neighbours(m,n)
    Bm=spdiags([ones(m,1) , zeros(m,1) , ones(m,1)], [-1 0 1], m, m);
    Bn=spdiags([ones(n,1) , zeros(n,1) , ones(n,1)], [-1 0 1], n, n);
    out=kron(speye(n),Bm)+kron(Bn,speye(m));
end
