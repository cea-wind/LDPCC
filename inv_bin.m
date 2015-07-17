function [ out ] = inv_bin( in )
%INV_BIN Summary of this function goes here
%   Detailed explanation goes here
%   input: A sparse binary matrix;
%	output: The inverse matrix of input matrix;  
[m,n] = size(in);
if(m~=n)
   fprintf('m~=n\n');
   return ;
end
E = eye(m);
%%

for i = 1:m
    noneZerosIndex = find(in(:,i)); 
    noneZerosIndex = noneZerosIndex(find(noneZerosIndex>=i));
    if(length(noneZerosIndex)==0)  
        randIndex = randi([i+1,m],1);
        temp = in(:,i);
        in(:,i) = in(:,randIndex);
        in(:,randIndex) = temp;       
        temp = E(:,i);
        E(:,i) = E(:,randIndex);
        E(:,randIndex) = temp;
    end
    id1 = noneZerosIndex(1);
    temp = in(i,:);
    in(i,:) = in(id1,:);
    in(id1,:) = temp;
    temp = E(i,:);
    E(i,:) = E(id1,:);
    E(id1,:) = temp;
    
    noneZerosIndex = find(in(:,i)); 
    for cc = 1:length(noneZerosIndex)
        if(noneZerosIndex(cc)~=i)  
            temp = mod(in(noneZerosIndex(cc),:)+in(i,:) , 2);
            in(noneZerosIndex(cc),:) = temp;
            temp = mod(E(noneZerosIndex(cc),:)+E(i,:) , 2);
            E(noneZerosIndex(cc),:) = temp;
        end
    end
end
out = E;
end

