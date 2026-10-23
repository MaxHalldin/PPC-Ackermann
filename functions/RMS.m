function [res] = RMS(X)
%RMS calcualtes RMS of X
    %X is array of scalars or vectors
arguments (Input)
    X
end

arguments (Output)
    res
end

temp = zeros(size(X,1),1);
for i=1:size(X,1)
    temp(i) = X(:,i).'*X(:,i);
end
res = sqrt(1/length(X) * sum(temp) );

end