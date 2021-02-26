% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
% Approximates tight extension with the procedure described in the paper
% (by increasing p)
% INPUT:
% g         - adjacency matrix of graph (size NxN, N number of vertices)
% f         - values on given vertices
% maske     - area on which the extension should be computed
% p         - maximal p of the Laplacians (capped at 2400)
% OUTPUT:
% f         - extended function values
% Llex      - Llex Vector
function [f,Llex]=ApproximateTight(g,f,maske,p_max)
% Set parameters
eps = 1e-10;
maxIter = 1000;

% Restrict graph to relevant part
relevant=((g*maske+maske)~=0);
graph_relevant=(g(relevant,relevant));
f_n=f(relevant,:);
maske_relevant=maske(relevant);
graph2 = logical(graph_relevant);
weights = graph_relevant;
weights = full(weights(weights>0).^2);
weights = weights/max(weights);
clearvars graph_relevant;

% Build difference matrix
[n,dim]=size(f_n);
k=1;
[ii,jj] = find(graph2);
nn = length(ii);
indices = zeros(2*nn,2);
values = zeros(2*nn,1);
for i=1:nn
    indices(k,:) = [i, ii(i)];
    values(k) = 1;
    k = k+1;
    indices(k,:) = [i, jj(i)];
    values(k) = -1;
    k = k+1;
end
A = sparse(indices(:,1),indices(:,2),values,nn,n);
clearvars ii jj values;

% Build export matrix
tmp = sum(graph2);
[num_edges, ~] = size(A);
inner_vert = nnz(maske_relevant);
edge = 1;
vert = 1;
count = 1;
indices = zeros(num_edges,2);
for i =1:n
    if maske_relevant(i)
        indices(count:count + tmp(i)-1,1)=vert;
        indices(count:count + tmp(i)-1,2)= edge:edge + tmp(i)-1;
        edge = edge + tmp(i);
        count = count + tmp(i);
        vert= vert+1;
    else
        edge = edge + tmp(i);
    end
end
indices = indices(1:count-1,:);
export = sparse(indices(:,1),indices(:,2),ones(count-1,1),...
    inner_vert,num_edges);
clearvars tmp indices;

% Compute logical vector with inner edges
[ind2, ind3] = find(graph2(maske_relevant,maske_relevant));
ind1 = graph2(:);
graph2(~maske_relevant,:) = false;
graph2(:,~maske_relevant) = false;
graph2 = graph2(:);
ind1 = full(graph2(ind1));
num_inner = sum(ind1);
clearvars graph2;

% Set index vectors for Hessian computation
indi2 = zeros(dim^2*(num_inner+inner_vert),1);
indi3 = zeros(dim^2*(num_inner+inner_vert),1);
count = 1;
part2 = dim^2*num_inner;
for i=1:dim
    for j=1:dim
        indi2(1+(count-1)*num_inner:count*num_inner,1) = ind2 +...
            (j-1)*inner_vert;
        indi3(1+(count-1)*num_inner:count*num_inner,1) = ind3 +...
            (i-1)*inner_vert;
        indi2(1+part2+(count-1)*inner_vert:part2+count*inner_vert,1) = ...
            (1:inner_vert)' + (j-1)*inner_vert;
        indi3(1+part2+(count-1)*inner_vert:part2+count*inner_vert,1) = ...
            (1:inner_vert)' + (i-1)*inner_vert;
        count = count + 1;
    end
end
indi2 = cast(indi2,'uint32');
indi3 = cast(indi3,'uint32');
[~,idx] = sortrows([indi2,indi3]);
idx = cast(idx,'uint32');
% Sort the indices
indi2 = indi2(idx);
indi3 = indi3(idx);
clearvars ind2 ind3;

% Solve problem for p=2
Llex = zeros(inner_vert,200);
counter= 1;
f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,indi3,idx,...
    maske_relevant,2,1,eps,maxIter);
f(relevant,:) = f_n;
Llex(:,counter) = computeL_lex(g,~maske,f);

% Multiscale procedure until maximal p is reached, capped at p = 2400
for p=5:5:20
    if p>p_max
        break;
    end
    disp(p);
    f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,...
        indi3,idx,maske_relevant,p,0.9,eps,maxIter);
    counter = counter +1;
    f(relevant,:) = f_n;
    Llex(:,counter) = computeL_lex(g,~maske,f);
end
for p=30:10:100
    if p>p_max
        break;
    end
    disp(p);
    f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,...
        indi3,idx,maske_relevant,p,0.95,eps,maxIter);
    counter = counter +1;
    f(relevant,:) = f_n;
    Llex(:,counter) = computeL_lex(g,~maske,f);
end
for p=  110:10:200
    if p>p_max
        break;
    end
    disp(p);
    f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,...
        indi3,idx,maske_relevant,p,0.978,eps,maxIter);
    counter = counter +1;
    f(relevant,:) = f_n;
    Llex(:,counter) = computeL_lex(g,~maske,f);
end
for p=210:10:400
    if p>p_max
        break;
    end
    disp(p);
    f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,...
        indi3,idx,maske_relevant,p,0.985,eps,maxIter);
    counter = counter +1;
    f(relevant,:) = f_n;
    Llex(:,counter) = computeL_lex(g,~maske,f);
end
for p=410:10:700
    if p>p_max
        break;
    end
    disp(p);
    f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,...
        indi3,idx,maske_relevant,p,0.99,eps,maxIter);
    counter = counter +1;
    f(relevant,:) = f_n;
    Llex(:,counter) = computeL_lex(g,~maske,f);
end
for p=710:10:1000
    if p>p_max
        break;
    end
    disp(p);
    f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,...
        indi3,idx,maske_relevant,p,0.9945,eps,maxIter);
    counter = counter +1;
    f(relevant,:) = f_n;
    Llex(:,counter) = computeL_lex(g,~maske,f);
end
for p=1010:10:1400
    if p>p_max
        break;
    end
    disp(p);
    f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,...
        indi3,idx,maske_relevant,p,0.9975,eps,maxIter);
    counter = counter +1;
    f(relevant,:) = f_n;
    Llex(:,counter) = computeL_lex(g,~maske,f);
end
for p=1410:10:1900
    if p>p_max
        break;
    end
    disp(p);
    f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,...
        indi3,idx,maske_relevant,p,0.9985,eps,maxIter);
    counter = counter +1;
    f(relevant,:) = f_n;
end
for p=1910:10:2400
    if p>p_max
        break;
    end
    disp(p);
    f_n=iterations_descent(f_n,A,export,weights,ind1,indi2,...
        indi3,idx,maske_relevant,p,0.9995,eps,maxIter);
    counter = counter +1;
    f(relevant,:) = f_n;
end
f(relevant,:) = f_n;
end

% Routine for minimization of E_p
function out = iterations_descent(f_n,A,export,weights,ind1,indi2,...
    indi3,idx,maske_relevant,p,tol,eps,maxIter)
%Initialize
GPU = false;%logical(gpuDeviceCount);
[~,dim]=size(f_n);
[num_edges, ~] = size(A);
inner_vert = nnz(maske_relevant);

% Compute solution explicity if p=2
if p==2
    out = f_n;
    LapOp = export*A(:,maske_relevant)+1e-7*speye(inner_vert);
    rhs = - export*A*(f_n.*(~maske_relevant));
    [L,U] = ilu(LapOp);
    for i=1:dim
        out(maske_relevant,i) = bicgstab(LapOp,rhs(:,i),1e-8,100,L,U);
    end
    return;
end

% Compute first Laplacian and build rescaling weight
das_im_p_lap = zeros(num_edges,dim);
f_diff=A*f_n;
f_diff_norm_sq = weights.*sum(f_diff.^2,2);
% Build sacle
scale = zeros(num_edges,1);
export2 = logical(export');
for i =1:inner_vert
    scale(export2(:,i),1) = ...
        max(f_diff_norm_sq(export2(:,i)));
end
scale(scale<1e-12) = 1;
scale = tol^2./scale;
% Compute factors
fac = powq(scale.*f_diff_norm_sq,(p-4)/2).*weights;
first_fac=fac.*scale.*f_diff_norm_sq;
for i=1:dim
    das_im_p_lap(:,i)=first_fac.*f_diff(:,i);
end
Lap=reshape(export*das_im_p_lap,inner_vert,dim);
clearvars f_diff_norm_sq das_im_p_lap;
infLap_norm=sum(sqrt(sum(Lap.^2,2)));
iterations=0;
abbruchkriterium=infLap_norm/inner_vert;

% Precompute some things
if GPU
    indP = zeros(dim,1);
    for j=1:dim
        indP(j,1) = 1 + (j-1)*(dim+1);
    end
end

% Iterate until stopping criterion is satisfied
while abbruchkriterium>eps && iterations<maxIter
    % Gradient direction
    descent = Lap;
    
    % Newton direction
    if p>4
        % Compute the Entries of the Hessian blockwise
        fac = (p-2)*scale.*fac.*weights;
        count = 1;
        Entries = zeros(num_edges,dim^2);
        for i = 1:dim
            for j = 1:dim
                if i<=j
                    if i == j
                        Entries(:,count) = first_fac +...
                            fac.*f_diff(:,i).*f_diff(:,j);
                        count = count + 1;
                    else
                        Entries(:,count) = fac.*f_diff(:,i).*f_diff(:,j);
                        count = count + 1;
                    end
                else
                    Entries(:,count) = Entries(:, (j-1)*dim + i);
                    count = count + 1;
                end
            end
        end
        clearvars f_diff fac first_fac;
        % Compute the diagonals of the blocks and adapt entries
        dia_ent = export*Entries;
        tmp = [-reshape(Entries(ind1,:),[],1); dia_ent(:)];
        clearvars Entries;
        % Delete small entries
        tmp(abs(tmp)<1e-9) = 0;
        % Build and solve Hessian
        if ~GPU
            % Cpu Code
            clearvars dia_ent;
            Hess = sparse(cast(indi3,'double'),cast(indi2,'double'),...
                tmp(idx),dim*inner_vert,dim*inner_vert);
            % Make Hessian positive definite
            Hess = Hess + 1e-7*speye(dim*inner_vert);
            D = tril(Hess);
            descent = reshape(...
                bicgstab(Hess,descent(:),1e-7,1000,D),inner_vert,dim);
        else
            %GPU Code
            % Build Hessian
            Hess = sparse(cast(indi3,'double'),cast(indi2,'double'),...
                gpuArray(tmp(idx)),dim*inner_vert,dim*inner_vert);
            Hess = Hess + 1e-7*speye(dim*inner_vert);
            % Build preconditioner
            tmp = dia_ent(:,indP);
            clearvars dia_ent
            tmp = tmp(:) + 1e-7* ones(dim*inner_vert,1);
            tmp= gpuArray(spdiags(tmp.^(-1),0,dim*inner_vert,dim*inner_vert));
            D = @(x) tmp*x;
            descent = gpuArray(descent);
            descent = reshape(...
                bicgstab(Hess,descent(:),1e-7,1500,D),inner_vert,dim);
            descent = gather(descent);
        end
        clearvars tmp Hess D;
    end
    f_n(maske_relevant,:) = f_n(maske_relevant,:) + descent;
    
    % Recompute Laplace
    das_im_p_lap = zeros(num_edges,dim);
    f_diff=A*f_n;
    f_diff_norm_sq = sum(f_diff.^2,2).*weights;
    
    %Recompute preconditioner if necessary
    if max(scale.*f_diff_norm_sq)^p>5
        scale = zeros(num_edges,1);
        for i=1:inner_vert
            scale(export2(:,i),1) = ...
                max(f_diff_norm_sq(export2(:,i)));
        end
        scale(scale<1e-12) = 1;
        scale = tol^2./scale;
    end
    fac = powq(scale.*f_diff_norm_sq,(p-4)/2).*weights;
    first_fac=fac.*scale.*f_diff_norm_sq;
    for i=1:dim
        das_im_p_lap(:,i)=first_fac.*f_diff(:,i);
    end
    Lap=reshape(export*das_im_p_lap,inner_vert,dim);
    clearvars f_diff_norm_sq das_im_p_lap;
    
    % Check convergence criterion
    infLap_norm=sum(sqrt(sum(Lap.^2,2)));
    iterations=iterations+1;
    abbruchkriterium=infLap_norm/inner_vert;
    if mod(iterations,50)==0
        disp(abbruchkriterium);
        disp(iterations);
    end
end
disp(abbruchkriterium)
disp(iterations)
out = f_n;
end
function out = powq(a,p)
if floor(p)==p && (p<101)
    out = a;
    for i =1:p-1
        out = out.*a;
    end
elseif floor(p)==p-0.5 && (p<101)
    out = sqrt(a);
    for i =1:floor(p)
        out = out.*a;
    end
else
    out = a.^p;
end
end
