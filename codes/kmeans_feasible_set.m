function sets = kmeans_feasible_set(M,snum,kdim)
% k-means clustering for selecting a subset of chips from a given set
% Input
% M - matrix of faithful from Walsh Hadamard matrix
% score - a vector with the Walsh score of each row of M
% snum - cardinality of set (i.e. number of clusters => k=snum
% kdim - dimension of embedding
% Output
% sets - 

[~, rsize]=size(M); % number of rows of the matrix
if kdim > rsize
    disp('size of lower-dimension has to be less than the length of a row of the input matrix')
    return;
end

lind=(1:1:rsize)';
M = M*randn(rsize,kdim)/sqrt(kdim); % "lower-dimensional" Gaussian embedding
%opts=statset('Display','final');
% Define the number of clusters and perform k-means clustering
%[idx, centroids] = kmeans(M, snum,Distance="sqeuclidean",Replicates=kdim,Options=opts);
[idx, centroids] = kmeans(M, snum,Distance="sqeuclidean",Replicates=kdim);
% Extract the subset of chips based on cluster indices
sets = zeros(snum,1);
for i = 1:snum
    ii = lind(idx==i);
    % In cluster i, find the element closest to the centroid
    [~,mm]=min(vecnorm(M(ii,:)-centroids(i,:),2,2)); 
    sets(i)=ii(mm);
end




