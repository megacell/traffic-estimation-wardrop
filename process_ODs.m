clear all; clc

load ODs/grouped_ODs.csv
A = sortrows(grouped_ODs);

row = 1;
i = 1;
while row<size(A,1)
    B(i,:) = A(row,:);
    k = 0;
    while sum(A(row,1:2)==A(row+k+1,1:2))>=2
        B(i,3) = B(i,3)+A(row+k+1,3);
        k = k+1;
    end
    i = i+1;
    row = row+k+1;
end
B(i,:) = A(row,:);

csvwrite('ODs/processed_ODs.csv',B)