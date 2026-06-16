function W = HadtoW(n)
% Hadamard matrices : - conversion from natural ordering to 
% sequency ordering (i.e. Walsh Hadamard matrix)
% Input
% n: order of matrix
% Output
% W: Walsh Hadamard matrix of order n
% Algorithm adapted from:
% Chen, C. F. and Leung, W. K., Algorithms for converting sequency-,
% dyadic-, and Hadamard-ordered Walsh functions, Mathematics and Computers
% in Simulation 27, 471-478, (1985)

% Generate index vector
twon=2^n;
id = zeros(twon,1);
for i=1:n
    a = 2^(n-i);
    b = 2^(i-1) - 0.5;
    for j = 0.5:1:b
        id(round(b+j)+1) = id(round(b-j)+1)+a;
    end
end
id=id+1;
H=hadamard(twon);
W=H(id,:);

