function[A,r,nul] = echelon(B) %% B is a matrix of any size. This function converts the matrix to row echelon form. Outputs a matrix A which is in row echelon form,r which is the rank of the matrix and nul is the nullity.
m = size(B,1); %% number of rows
n = size(B,2); %% number of columns
A = B;
c = zeros(1,n);
r = m;
for j = 1:m-1
    x = A(j,j);
    if x ~=0
        A(j,:) = A(j,:)./x;                                     
    end
    for i = j+1:m
        y = A(i,j);                                               
        A(i,:) = A(i,:) - (y*A(j,:));
    end
end
for i = 1:m-1
    for j = i+1:m
        if((isequal(A(i,:),c)==1) & ((isequal(A(j,:),c)==0)))
            temp = A(j,:);
            temp2 = A(i,:);
            A(i,:) = temp;
            A(j,:) = temp2;                  
        end
    end
end
for i = 1:m
    for j = 1:n
        if(mod(A(i,j),1) ~= 0)
            A(i,:) = A(i,:).*120;            %% To make sure that no decimals are present
            break
        end
    end
end
for i = 1:m
    x = max(abs(A(i,:)));
    for j = 1:n
        y = A(i,j);                          %% Simplifies matrix
        if(y>0)
            x = gcd(x,y);
        end
        
    end
    if(x>0)
        A(i,:) = A(i,:)./x;
    end
end
for i = 1:m-1
    for j = i+1:m
        d = A(i,:) - A(j,:);
        if((isequal(d,c)==1))
            A(j,:) = c;
        end
    end
end
    
for i = 1:m
    if((isequal(A(i,:),c)==1))              %% Finds rank of matrix
        r = r-1;
    end
end
nul = n-r;
end
        