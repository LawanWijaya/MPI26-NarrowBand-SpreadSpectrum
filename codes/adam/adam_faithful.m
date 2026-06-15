function return_data = adam_faithful(inp_seq,cutoff)

    S = [0; cumsum(inp_seq)];

    M = length(S);
    k = [0:floor(M/2)-1 -ceil(M/2):-1]';
    xi = 2*pi*k;

    Fhat = fft(S);
    mask = abs(xi) <= cutoff;
    f_filt = real(ifft(Fhat .* mask));

    seq_new = 2*(diff(f_filt) >= 0) - 1;

    %one = figure(1);
    %clf(one);
    %one_ax = axes(one);
    %plot(one_ax,1:length(S),S);
    %title(one_ax,'Input');

    %two = figure(2);
    %clf(two);
    %two_ax = axes(two);
    %plot(two_ax,1:length(f_filt),f_filt);
    %title(two_ax,'Output');

    return_data = seq_new;

end