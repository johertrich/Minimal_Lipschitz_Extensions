% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
function out=refine_solution(graph,U,erg,n)
% The algorithm tries to optimize the Lipschitz constant for every vertex
% in a random fashion. The procedure is meant as a second step, since it
% does not necesarilly converge to the tight solution. However, the result
% is always tighter as the input.
% The scheme can be found in Section 4 of the PhD Thsis 'Minimall Lipschitz
% extensions' by Phan Thanh.
% INPUT:
%   graph - adjacency matrix of the graph (size NxN, where N is the number
%           of the vertices)
%   U     - logical if vertex i is in the boundary (size Nx1)
%   erg   - vertex values produced by ADMM(size N*nx1, where n is the
%           dimension of each value at one vertex)
%   n     - dimension of the problem
% OUTPUT:
%   out   - refined soltion

% Intialize output;
out=erg;
% Reorder input
erg2 = zeros(length(erg)/n,n);
for i=1:n
    erg2(:,i) = erg(i:n:end);
end
% Determine non boundary vertices
free = 1:length(U);
free = free(~U);
% Random order for visiting vertices
per = randperm(length(free));
S_max = zeros(1,length(free));
for m=1:length(free)
    % Break paramater for exiting after the optimal vertex value is
    % determined
    b_para=false;
    i = free(per(m));
    % Adjacent vertex values
    adj_ind = (graph(:,i)~= 0);
    adj_weights = 1./full(graph(adj_ind,i));
    tmp = erg2(adj_ind,:);
    len = size(tmp,1);
    % Old value of S
    S_old = max(sum((tmp-repmat(erg2(i,:),len,1)).^2,2)./(adj_weights.^2));
    S_max(i)=S_old;
    %
    % First stage of procedure (convex combination 2 points)
    %
    tmp2 = reshape(tmp,[1,len,n]);
    a1 = repmat(tmp2,len,1,1);
    w2 = repmat(adj_weights',len,1,1);
    a2 = permute(a1,[2,1,3]);
    w1 = permute(w2,[2,1,3]);
    % Optimal vertex value between vertex i and j
    fij = (w1.*a1 +w2.*a2)./(w1 + w2);
    % Optimal Lipschitz bound for vertex i and j
    Lij = sum((a1-a2).^2,3)./((w1(:,:,1) + w2(:,:,1)).^2);
    % Check if the Lipschitz bound together with the optimal vertex value
    % satisfies S(free(per(m))) \leq Lij
    test = zeros(size(Lij));
    for j= 1:len
        bound = sum((fij - repmat(tmp2(1,j,:),len,len,1)).^2,3)<=...
            (Lij*adj_weights(j)^2+1e-15);
        bound(:,j) = true;
        bound(j,:) = true;
        test = test + bound;
    end
    % If more values are possible, make a random choice
    test = triu(test,1);
    [ind1, ind2] = find(test==len);
    if length(ind1)>2
        ind = randi(length(ind1));
        ind1 = ind1(ind);
        ind2 = ind2(ind);
    end
    % Write optimal solution into input and output
    if ~isempty(ind1)
        out(1 +(i-1)*n:i*n) = fij(ind1(1),ind2(1),:);
        S_new = max(sum((tmp-repmat(out(1 +(i-1)*n:i*n)',len,1)).^2,...
            2)./(adj_weights.^2));
        if S_old-S_new>1e-12
            erg2(i,:) =  out(1 +(i-1)*n:i*n);
            % Go to next vertex
            b_para = true;
        else
            out(1 +(i-1)*n:i*n) = erg2(i,:);
        end
    end
    %
    % Second stage of procedure (convex combination 3 points)
    %
    if ~b_para&&(len>2)&&(n>1)
        % Check all triplets of vertices in a random order
        indices = nchoosek(1:len,3);
        s_ind = size(indices,1);
        r_ind = randperm(s_ind);
        for ind = 1:s_ind
            j= indices(r_ind(ind),1);
            k= indices(r_ind(ind),2);
            l= indices(r_ind(ind),3);
            % Compute determinant, weights are added
            ajk = sum((tmp(j,:)-tmp(k,:)).^2);
            ajl = sum((tmp(j,:)-tmp(l,:)).^2);
            akl = sum((tmp(k,:)-tmp(l,:)).^2);
            wj = adj_weights(j)^2;
            wk = adj_weights(k)^2;
            wl = adj_weights(l)^2;
            cc = -2*ajk*ajl*akl;
            bb = 2*ajk*ajl*wk - 2*akl^2*wj - 2*ajk^2*wl - 2*ajl^2*wk +...
                2*ajk*akl*wj + 2*ajl*akl*wj + 2*ajk*ajl*wl +...
                2*ajl*akl*wk + 2*ajk*akl*wl;
            aa = 2*ajl*wj*wk - 2*akl*wj^2 - 2*ajk*wl^2 - 2*ajk*wj*wk -...
                2*ajl*wk^2 + 2*ajk*wj*wl - 2*ajl*wj*wl + 2*akl*wj*wk +...
                2*ajk*wk*wl + 2*ajl*wk*wl + 2*akl*wj*wl - 2*akl*wk*wl;
            % Compute solution of polynomial, which is the
            % optimal Lipschitz constant
            sol = solve_eq(aa,bb,cc);
            % Check determinant condition and solvability
            if ~isempty(sol)
                for ii=1:length(sol)
                    % Compute intersection of 3 spheres
                    [p_12_a,~] = trilaterate(tmp(j,:),tmp(k,:),tmp(l,:),...
                        sol(ii)*wj,sol(ii)*wk,sol(ii)*wl,[],[]);
                    % Check solution
                    if ~isempty(p_12_a)
                        % Shorten solution if n=2
                        p_12_a = p_12_a(1:n);
                        % Check if S(free(per(m))) \leq sol
                        check = sum((repmat(p_12_a,len,1)-tmp).^2,2)./(...
                            adj_weights.^2)<=(sol(ii)+1e-15);
                        check([j,k,l]) = true;
                        check = sum(check);
                        if check==len
                            % Check if p_12_a is in the convex hull was done in
                            % trilaterate. Compute new value of S and double-
                            % check with old value
                            S_new = max(sum((tmp-repmat(...
                                p_12_a,len,1)).^2,2)./(adj_weights.^2));
                            if S_old-S_new>1e-12
                                % Write optimal solution into input
                                % and output
                                out(1 +(i-1)*n:i*n) = p_12_a;
                                erg2(i,:) = p_12_a;
                                % Go to next vertex
                                b_para = true;
                                break;
                            end
                        end
                    end
                end
                if b_para
                    break;
                end
            end
        end
    end
    %
    % Third stage of procedure (convex combination of 4 points)
    %
    if ~b_para&&(n>2)&&(len>3)
        % Check all quadruples of vertices in random order
        indices = nchoosek(1:len,4);
        s_ind = size(indices,1);
        r_ind = randperm(s_ind);
        for ind = 1:s_ind
            j= indices(r_ind(ind),1);
            k= indices(r_ind(ind),2);
            l= indices(r_ind(ind),3);
            o= indices(r_ind(ind),4);
            % Compute determinant, weights are added
            ajk = sum((tmp(j,:)-tmp(k,:)).^2);
            ajl = sum((tmp(j,:)-tmp(l,:)).^2);
            ajo = sum((tmp(j,:)-tmp(o,:)).^2);
            akl = sum((tmp(k,:)-tmp(l,:)).^2);
            ako = sum((tmp(k,:)-tmp(o,:)).^2);
            alo = sum((tmp(l,:)-tmp(o,:)).^2);
            wj = adj_weights(j)^2;
            wk = adj_weights(k)^2;
            wl = adj_weights(l)^2;
            wo = adj_weights(o)^2;
            cc = - ajk^2*alo^2 + 2*ajk*ajl*ako*alo + 2*ajk*ajo*akl*alo -...
                ajl^2*ako^2 + 2*ajl*ajo*akl*ako - ajo^2*akl^2;
            bb = 2*ajl*ako^2*wj + 2*ajo*akl^2*wj + 2*ajk*alo^2*wj +...
                2*ajl^2*ako*wk + 2*ajo^2*akl*wk + 2*ajk*alo^2*wk +...
                2*ajl*ako^2*wl + 2*ajo^2*akl*wl + 2*ajk^2*alo*wl +...
                2*ajo*akl^2*wo + 2*ajl^2*ako*wo + 2*ajk^2*alo*wo -...
                2*ajl*ajo*akl*wk - 2*ajl*akl*ako*wj - 2*ajl*ajo*ako*wk -...
                2*ajo*akl*ako*wj - 2*ajk*ajl*ako*wl - 2*ajk*ajl*alo*wk -...
                2*ajk*ajo*akl*wl - 2*ajk*akl*alo*wj + 4*ajk*ajo*ako*wl -...
                2*ajk*ajo*alo*wk - 2*ajk*ako*alo*wj - 2*ajl*ajo*ako*wl +...
                4*ajl*ajo*alo*wk - 2*ajl*ako*alo*wj - 2*ajo*akl*alo*wj -...
                2*ajk*ajo*alo*wl - 2*ajl*ako*alo*wk - 2*ajo*akl*ako*wl -...
                2*ajo*akl*alo*wk + 4*akl*ako*alo*wj - 2*ajk*ako*alo*wl +...
                4*ajk*ajl*akl*wo - 2*ajk*ajl*ako*wo - 2*ajk*ajo*akl*wo -...
                2*ajl*ajo*akl*wo - 2*ajk*ajl*alo*wo - 2*ajl*akl*ako*wo -...
                2*ajk*akl*alo*wo;
            aa = - ajk^2*wl^2 + 2*ajk^2*wl*wo - ajk^2*wo^2 + 2*ajk*ajl*wk*wl -...
                2*ajk*ajl*wk*wo - 2*ajk*ajl*wl*wo + 2*ajk*ajl*wo^2 -...
                2*ajk*ajo*wk*wl + 2*ajk*ajo*wk*wo + 2*ajk*ajo*wl^2 -...
                2*ajk*ajo*wl*wo + 2*ajk*akl*wj*wl - 2*ajk*akl*wj*wo -...
                2*ajk*akl*wl*wo + 2*ajk*akl*wo^2 - 2*ajk*ako*wj*wl +...
                2*ajk*ako*wj*wo + 2*ajk*ako*wl^2 - 2*ajk*ako*wl*wo +...
                4*ajk*alo*wj*wk - 2*ajk*alo*wj*wl - 2*ajk*alo*wj*wo -...
                2*ajk*alo*wk*wl - 2*ajk*alo*wk*wo + 4*ajk*alo*wl*wo -...
                ajl^2*wk^2 + 2*ajl^2*wk*wo - ajl^2*wo^2 + 2*ajl*ajo*wk^2 -...
                2*ajl*ajo*wk*wl - 2*ajl*ajo*wk*wo + 2*ajl*ajo*wl*wo +...
                2*ajl*akl*wj*wk - 2*ajl*akl*wj*wo - 2*ajl*akl*wk*wo +...
                2*ajl*akl*wo^2 - 2*ajl*ako*wj*wk + 4*ajl*ako*wj*wl -...
                2*ajl*ako*wj*wo - 2*ajl*ako*wk*wl + 4*ajl*ako*wk*wo -...
                2*ajl*ako*wl*wo - 2*ajl*alo*wj*wk + 2*ajl*alo*wj*wo +...
                2*ajl*alo*wk^2 - 2*ajl*alo*wk*wo - ajo^2*wk^2 +...
                2*ajo^2*wk*wl - ajo^2*wl^2 - 2*ajo*akl*wj*wk - 2*ajo*akl*wj*wl +...
                4*ajo*akl*wj*wo + 4*ajo*akl*wk*wl - 2*ajo*akl*wk*wo -...
                2*ajo*akl*wl*wo + 2*ajo*ako*wj*wk - 2*ajo*ako*wj*wl -...
                2*ajo*ako*wk*wl + 2*ajo*ako*wl^2 - 2*ajo*alo*wj*wk +...
                2*ajo*alo*wj*wl + 2*ajo*alo*wk^2 - 2*ajo*alo*wk*wl - akl^2*wj^2 +...
                2*akl^2*wj*wo - akl^2*wo^2 + 2*akl*ako*wj^2 - 2*akl*ako*wj*wl -...
                2*akl*ako*wj*wo + 2*akl*ako*wl*wo + 2*akl*alo*wj^2 -...
                2*akl*alo*wj*wk - 2*akl*alo*wj*wo + 2*akl*alo*wk*wo -...
                ako^2*wj^2 + 2*ako^2*wj*wl - ako^2*wl^2 + 2*ako*alo*wj^2 -...
                2*ako*alo*wj*wk - 2*ako*alo*wj*wl + 2*ako*alo*wk*wl -...
                alo^2*wj^2 + 2*alo^2*wj*wk - alo^2*wk^2;
            % Compute solution of polynomial, which is the
            % optimal Lipschitz constant
            sol = solve_eq(aa,bb,cc);
            % Check simplex condtion
            safe = abs(det([tmp(k,:)-tmp(j,:);tmp(l,:)-tmp(j,:);...
                tmp(o,:)-tmp(j,:)]))>1e-14;
            if ~isempty(sol) && safe
                for ii=1:length(sol)
                    % Compute intersection of 4 spheres
                    [p_12_a,p_12_b] = trilaterate(tmp(j,:),tmp(k,:),...
                        tmp(l,:),sol(ii)*wj,sol(ii)*wk,sol(ii)*wl,...
                        tmp(o,:),sol(ii)*wo);
                    % Check first solution
                    if ~isempty(p_12_a)
                        % Check if S(free(per(m))) \leq sol
                        check = sum((repmat(p_12_a,len,1)-tmp).^2,2)./(...
                            adj_weights.^2)<=(sol(ii)+1e-15);
                        check([j,k,l,o]) = true;
                        check = sum(check);
                        if check==len
                            % Check if p_12_a is in the convex hull
                            Mat = [tmp(k,:)-tmp(j,:);tmp(l,:)-tmp(j,:);...
                                tmp(o,:)-tmp(j,:)]';
                            lhs = (p_12_a-tmp(j,:))';
                            paras = Mat\lhs;
                            check = (norm(Mat*paras-lhs) <1e-14);
                            in_hull = (sum(paras)<(1+1e-15))&...
                                (sum(paras>-1e-15)==3)&check;
                            if in_hull
                                % Compute new value of S and double-
                                % check with old value,
                                S_new = max(sum((tmp-repmat(...
                                    p_12_a,len,1)).^2,2)./(adj_weights.^2));
                                if S_old-S_new>1e-12
                                    % Write optimal solution into input
                                    % and output
                                    out(1 +(i-1)*n:i*n) = p_12_a;
                                    erg2(i,:) = p_12_a;
                                    % Go to next vertex
                                    b_para= true;
                                    break;
                                end
                            end
                        end
                    end
                    % Check second solution, similar as previous
                    % case.
                    if ~isempty(p_12_b)
                        check = sum((repmat(p_12_b,len,1)-tmp).^2,2)./(...
                            adj_weights.^2)<=(sol(ii)+1e-15);
                        check([j,k,l,o]) = true;
                        check = sum(check);
                        if check==len
                            Mat = [tmp(k,:)-tmp(j,:);tmp(l,:)-tmp(j,:);...
                                tmp(o,:)-tmp(j,:)]';
                            lhs = (p_12_b-tmp(j,:))';
                            paras = Mat\lhs;
                            check = (norm(Mat*paras-lhs) <1e-14);
                            in_hull = (sum(paras)<(1+1e-15))&...
                                (sum(paras>-1e-15)==3)&check;
                            if in_hull
                                S_new = max(sum((tmp-repmat(...
                                    p_12_b,len,1)).^2,2)./(adj_weights.^2));
                                if S_old-S_new>1e-12
                                    out(1 +(i-1)*n:i*n) = p_12_b;
                                    erg2(i,:) = p_12_b;
                                    b_para = true;
                                    break;
                                end
                            end
                        end
                    end
                end
                if b_para
                    break;
                end
            end
        end
    end
end
S_max = sort(S_max);
S_max((end-5):end);
end

function roots=solve_eq(a,b,c)
% Solves the quadratic equation a*x^2 +b*x +c = 0
% Case a \neq 0
if abs(a)>1e-15
    Disk = b^2-4*a*c;
    if Disk>0
        sol1 = (-b - sqrt(Disk))/(2*a);
        sol2 = (-b + sqrt(Disk))/(2*a);
        roots = [sol1 sol2];
    else
        roots = [];
    end
    % Case a = 0
else
    if abs(b)>1e-15
        sol1 = -c/b;
        roots = sol1;
    else
        roots = [];
    end
end
roots = roots(roots>=0);
end

function [p_12_a,p_12_b] = trilaterate(P1,P2,P3,r12,r22,r32,P4,r4)
% Computes the weighted Chebychev center of 3 (triangle) or 4 points (tetrahedron).
% Note that the squared radius is used as input. P4 and r4 can be empty.
if length(P1)==2
    P1 = [P1 0];
    P2 = [P2 0];
    P3 = [P3 0];
end
temp1 = P2-P1;
e_x = temp1/norm(temp1);
temp2 = P3-P1;
i = sum(e_x.*temp2);
temp3 = temp2 - i*e_x;
e_y = temp3/norm(temp3);
e_z = cross(e_x,e_y);
d = norm(P2-P1);
j = sum(e_y.*temp2);
x = (r12 - r22 + d*d)/(2*d);
y = (r12 - r32 -2*i*x + i*i + j*j)/(2*j);
temp4 = r12 - x*x - y*y;

% Here custom computations take place,
z = sqrt(abs(temp4));
if isempty(P4)
    p_12_a = P1 + x*e_x + y*e_y;
    % Check for to big values of z
    if z>5e-8
        p_12_a = [];
    end
    % Check convex combination
    coeffs = [(x - y*i/norm(temp3))/norm(temp1),y/norm(temp3)];
    if any(coeffs<-1e-15) || (sum(coeffs)>(1+1e-15))
        p_12_a = [];
    end
    % Second component not needed here
    p_12_b =[];
    return;
else
    % Check for imaginary root
    if temp4<0
        p_12_a=[];
        p_12_b=[];
        return;
    end
    % Write values
    p_12_a = P1 + x*e_x + y*e_y + z*e_z;
    p_12_b = P1 + x*e_x + y*e_y - z*e_z;
    % Check values
    if abs(norm(p_12_a-P4)^2-r4)>1e-13
        p_12_a =[];
    end
    if abs(norm(p_12_b-P4)^2-r4)>1e-13
        p_12_b =[];
    end
end
end
