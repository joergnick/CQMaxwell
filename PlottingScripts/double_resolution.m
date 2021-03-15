function [ Ad] = double_resolution( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n=length(A(1,:));
m=length(A(:,1));
Ad=zeros(2*m-1,2*n-1);
size(A)
Ad(1:2:2*m-1,1:2:2*n-1)=A;
Ad(1:2:2*m-1,2:2:2*n-2)=0.5*(Ad(1:2:2*m-1,1:2:2*n-3)+Ad(1:2:2*m-1,3:2:2*n-1));
Ad(2:2:2*m-2,:)=0.5*(Ad(1:2:2*m-3,:)+Ad(3:2:2*m-1,:));


end

