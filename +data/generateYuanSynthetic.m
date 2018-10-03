function [XU,S] = generateYuanSynthetic(p,h,r,n)
% construct data for synthetic data experiment in Chadrasekaran et al.
% 2012 - Ming Yuan Georgia Tech, section 6.1 - 36-cycle, 2 latent variables

%     p = number of observed variables
%     h = number of hidden variables
%     r = degree of connectivity
%     n = number of datapoints to generate

    S = eye(p+h);
    % p random locations on a 1x1 square
    locs = rand(p, 2);
    % Add partial correlation edges 
    for ii = 1:p
        for jj = 1:p
            if ii < jj
                add_edge = rand() < (2 * normpdf(norm(locs(ii,:) - locs(jj,:),2) * sqrt(p)));
                S(ii, jj) = r * add_edge;
            end
        end
    end
    % Remove edges when nodes have more than four edges connected to it
    max_edges = 2;
    for ii = 1:p
        [node_edges] = find(S(ii,:) > 0);
        if length(node_edges) > max_edges
            idx_delete = randperm(length(node_edges));
            S(ii, node_edges(idx_delete(max_edges+1:end))) = 0;
            S(node_edges(idx_delete(max_edges+1:end)),ii) = 0;
        end 
        S(ii, ii) = 1;
    end
    S = (S + S') - (eye(p+h).*diag(S)); 

    % Latent variables connected to random subset of observed variables
    subset_percent = 1;
    for ii = (p+1):(p+h)
        ind = 1:p;
        ind = randperm(p,round(subset_percent*p));
        vals = 0.11*rand(1,round(subset_percent*p));
        vals = 0.11*rand(1,round(subset_percent*p));
        S(ind,ii) = vals;
        S(ii, ind) = vals;
        for jj = (p+1):(p+h)
            if ii == jj
                S(ii, jj) = 1;
            else
                % Latent Variables connection weights
                S(ii, jj) = 0;
            end  
        end
    end

    C = inv(S);
    % XU is [p+h x n] synthetic data
    XU = (randn(n, p+h) * chol(C))';
end