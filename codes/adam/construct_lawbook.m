function return_data = construct_lawbook(w,faithful,find_best)

    cds = floor(length(faithful)/(2^w));

    lawbook = zeros(cds,w^2);

    for cd = 1:cds

        best = find_best(faithful,w);
        lawbook(cd,:) = best;
        faithful = faithful(~ismember(faithful,best));

    end

    return_data = lawbook;

end