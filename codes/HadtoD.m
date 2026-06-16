function D = HadtoD(n)
% Hadamard matrices : - conversion from natural ordering to 
% dyadic ordering 
% Input
% n: order of matrix
% Output
% D: Hadamard matrix of order n with dyadic ordering
% Algorithm adapted from:
% Chen, C. F. and Leung, W. K., Algorithms for converting sequency-,
% dyadic-, and Hadamard-ordered Walsh functions, Mathematics and Computers
% in Simulation 27, 471-478, (1985)

% Generate index vector
twon=2^n;
id = zeros(twon,1);
for i=1:n
    a = 2^(i-1);
    b = 2^(n-i);
    for j = 0:1:a-1
        id(a+j+1) = b+id(j+1);
    end
end
id=id+1;
H=hadamard(twon);
D=H(id,:);

