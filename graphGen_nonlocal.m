% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
% Computes an adjacency matrix of a graph on an image based on patch similarity
% INPUT:
% img               - image
% mask              - unknown pixels
% OPTIONAL:
% patchSize         - A patch around a pixel (i,j) is defined by 
%                     [i-patchSize:i+patchSize,j-patchSize:j+patchSize]
% radius            - Search similar patches in a neighborhood of size 2*radius+1
% needed_area       - Pixels, where the graph is generated
%                     Logical-Indexing-Matrix
% nNeighbors        - number of neighbors in the graph
% sigma             - weighting parameter. For smaller sigma there are more
%                     significant differences in the weights
% OUTPUT:
% out               - adjacency matrix of the resulting graph
function out=graphGen_nonlocal(img,mask,varargin)
    ip=inputParser;
    addOptional(ip,'patchSize',7);
    addOptional(ip,'radius',15);
    addOptional(ip,'needed_area',mask);
    addOptional(ip,'nNeighbors',8);
    addOptional(ip,'sigma',1);
    parse(ip,varargin{:});
    
    s=ip.Results;
    p=s.patchSize;
    r=s.radius;
    needed_area=s.needed_area;
    k=s.nNeighbors;
    sigma=s.sigma;
    out = graphGen_inner(img,mask,p,r,needed_area,k,sigma);
end

function out=graphGen_inner(img,mask,p,r,needed_area,k,sigma)
[m,n,d]=size(img);
img2=img.*repmat(1-mask,1,1,d); % Set unknown values to 0
% Insert zeros on the boundary
img2=[zeros(p+r,n+2*(p+r),d);zeros(m,p+r,d), img2,zeros(m,p+r,d);zeros(p+r,n+2*(p+r),d)];
% Expand mask
mask2=logical([ones(p+r,n+2*(p+r));ones(m,p+r),mask,ones(m,p+r);ones(p+r,n+2*(p+r))]);
[indx,indy]=find(needed_area==1);
nPixel=size(indx,1);
shift=p+r; % Shift because of the added boundary
out2=sparse(m*n,nPixel);
parfor i=1:nPixel
    
    if mod(i,500)==0
        disp('500 Pixels computed');
    end
    
    x=indx(i);
    y=indy(i);
    
    weights=zeros(2*r+1,2*r+1);
    nWorks=zeros(2*r+1,2*r+1);
    
    % All known points in patch around (x,y)
    [A,B]=find(mask2(x-p+shift:x+p+shift,y-p+shift:y+p+shift)==0);
    A=A-p-1;
    B=B-p-1;
    for ind=1:size(A,1)
        a=A(ind);
        b=B(ind);
        % Compare the pixels in the patches corresponding to (x,y) with the
        % pixels from the patches in radius r around (x,y)
        tmp_weight = double(~mask2(x+a-r+shift:x+a+r+shift,y+b-r+shift:y+b+r+shift));
        tmp = img2(x+a-r+shift:x+a+r+shift,y+b-r+shift:y+b+r+shift,:) -...
            repmat(img2(x+a+shift,y+b+shift,:),2*r+1,2*r+1,1);
        tmp = sum(tmp.^2,3);%.^0.5;
        weights=weights+tmp.*tmp_weight;
        nWorks=nWorks+tmp_weight; % Count comparable pixels
    end
    
    nWorks(nWorks==0) =1;
    weights=(1./nWorks).*weights;
    weights=exp(-weights/sigma^2).*(nWorks>=(2*p+1)^2/9); % Weighting
    weights(r+1,r+1)=0;
    largestWeights=zeros(2*r+1,2*r+1);
    for j=1:k % Compute the k largest weights
        max_value=max(max(weights));
        [x_max,y_max]=find(weights==max_value);
        len_max = length(x_max);
        ind2 = randi(len_max);
        x_max = x_max(ind2);
        y_max = y_max(ind2);
        largestWeights(x_max,y_max)=weights(x_max,y_max);
        weights(x_max,y_max)=0;
    end
    tempLine=sparse(m+2*r,n+2*r);
    tempLine(x:x+2*r,y:y+2*r)=largestWeights; % Build graph line for line
    tempLine=tempLine(r+1:m+r,r+1:n+r);
    out2(:,i)=reshape(tempLine,m*n,1);
end

out=sparse(m*n,m*n);
for i=1:nPixel
    x=indx(i);
    y=indy(i);
    out(:,x+m*(y-1))=sparse(out2(:,i));
end
out=max(out,out');
end
