function return_data = adam_faithful(inp_seq,cutoff)

    % choose code length
    code_length = length(inp_seq);

    % rudimentary discrete ADPCM (only at intervals)
    inp_con = [0;cumsum(inp_seq)];

    % fill in intervals
    m = 64;
    inp_con_fill = repelem(inp_con,m);
    M = length(inp_con_fill);

    % visualize cont input
    %inp_anal_fig = figure(1);
    %clf(inp_anal_fig);
    %inp_anal_ax = axes(inp_anal_fig);
    %plot(inp_anal_ax,1:M,inp_con_fill);
    %title(inp_anal_ax,'Input Analog');

    % Define frequencies and cutoff
    k = [0:M/2-1 -M/2:-1]';
    xi = 2*pi*k/code_length;

    % fft -> filter -> ifft
    Fhat = fft(inp_con_fill);
    mask = abs(xi) <= cutoff;
    Fhat_filt = Fhat .* mask;
    out_con_fill = real(ifft(Fhat_filt));

    % visualize cont output
    %out_anal_fig = figure(2);
    %clf(out_anal_fig);
    %out_anal_ax = axes(out_anal_fig);
    %plot(out_anal_ax,1:M,out_con_fill);
    %title(out_anal_ax,'Output Analog');

    % extract values at intervals
    out_cont_fill = out_con_fill(1:m:end);
    out_seq = 2*(diff(out_cont_fill) >= 0) - 1;

    return_data = out_seq;

    %inp_seq.'
    %out_seq.'

    %num_diff = sum(inp_seq ~= out_seq);

    %return_data = num_diff;

end