function return_data = adam_faithful(inp_seq,cutoff)
% ADAM_FAITHFUL applies a low-pass filter characterized by some given
% cutoff to a given input sequence, and returns the result
%
%last updated 06/26/2026 by Adam Petrucci
arguments
    inp_seq       % chip sequence to process
    cutoff        % truncation value for low-pass filter
end

    % Construct coarse-mesh difference representation of input sequence in
    % physical (as opposed to discrete) space. cumsum is built-in Matlab
    % function for constructing cumulative sum of vector
    S = [0; cumsum(inp_seq)];

    % Set up FFT
    M = length(S);
    k = [0:floor(M/2)-1 -ceil(M/2):-1]';
    xi = 2*pi*k;

    % Apply low-pass filter in frequency space using built-in Matlab
    % functionality for Fourier analysis
    Fhat = fft(S);
    mask = abs(xi) <= cutoff;
    f_filt = real(ifft(Fhat .* mask));

    % Interpret output in physical space as chip-sequence. diff computes
    % the difference between subsequent values
    seq_new = 2*(diff(f_filt) >= 0) - 1;

    % Return data
    return_data = seq_new;

end