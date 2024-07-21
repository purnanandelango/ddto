function costval = compute_ddto_cost(X,i)

    n = length(X);
    Ni = size(X{i},2);

    costval = 0;
    for j = setdiff(1:n,i)
        Nji = min(size(X{j},2),Ni);
        for k = 1:Nji
            % Percent difference
            if norm(X{j}(:,k) - X{i}(:,k))/norm(X{i}(:,k)) > 1e-3
                costval = costval + 1;
            end
        end
    end

end