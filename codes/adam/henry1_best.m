function return_data = henry1_best(faithful,w,options)
% HENRY1_BEST is a wrapper for applying Henry's algorithm for selecting an
% optimal codebook from a given set of faithful encodings.
%
%last updated 06/28/2026 by Adam Petrucci
arguments (Input)
    faithful                  % 
    w                         %
end
arguments (Input)
    options.kdim = 100        %
    options.n = 9             %
end

    Wn = HadtoW(options.n);
    M = Wn(faithful,:);
    snum = 2^w;
    kdim = options.kdim;

    henry_res = henry_best(M,faithful,snum,kdim,1);

    return_data = sort(faithful(henry_res));

end