% Draft
% Game of life
%n matrix dimension 
%spsize matrix sparse size
%tf final time

function tentativa_life(n,spsize,tf)
if nargin<3
    n=20;
    spsize=0.3;
    tf=100;
end

A=sprand(n,n,spsize);
B=spones(A);
%B=full(B);
% spy(B)
[row,col]=find(B);
pop=zeros(1,tf);
for t=1:tf
    pop(t)=sum(sum(B));
%Matriz B modification    
for i=1:length(row)

%% Treatment of edges

%First line
if row(i)==1 
    %First corner
    if col(i)==1
        neighbors = B(row(i)+1,col(i))+B(row(i),col(i)+1)+B(row(i)+1,col(i)+1);
        if neighbors ==3 || neighbors ==1
            B(row(i),col(i))=0;
        else
            B(row(i),col(i))=1;
        end
       
    %Second corner
    elseif col(i)==n
        neighbors = B(row(i)+1,col(i))+B(row(i),col(i)-1)+B(row(i)+1,col(i)-1);
        if neighbors ==3 || neighbors ==1
            B(row(i),col(i))=0;
        else
            B(row(i),col(i))=1;
        end
        
    else
        neighbors = B(row(i)+1,col(i))+B(row(i),col(i)-1)+B(row(i)+1,col(i)-1)+B(row(i)+1,col(i)+1)+B(row(i),col(i)+1);
        if neighbors >=3 || neighbors ==1
            B(row(i),col(i))=0;
        else
            B(row(i),col(i))=1;
        end
    end
    
%End line
elseif row(i)==n
    %Third corner
    if col(i)==1
        neighbors = B(row(i)-1,col(i))+B(row(i),col(i)+1)+B(row(i)-1,col(i)+1);
        if neighbors ==3 || neighbors ==1
            B(row(i),col(i))=0;
        else
            B(row(i),col(i))=1;
        end
        
        %Fourth corner
    elseif  col(i)==n
        neighbors = B(row(i)-1,col(i))+B(row(i),col(i)-1)+B(row(i)-1,col(i)-1);
        if neighbors ==3 || neighbors ==1
            B(row(i),col(i))=0;
        else
            B(row(i),col(i))=1;
        end
    else
        neighbors = B(row(i)-1,col(i))+B(row(i),col(i)-1)+B(row(i)-1,col(i)-1)+B(row(i)-1,col(i)+1)+B(row(i),col(i)+1);
        if neighbors >=3 || neighbors ==1
            B(row(i),col(i))=0;
        else
            B(row(i),col(i))=1;
        end
    end

%First column
elseif col(i)==1
    if row(i)~=1 && row(i)~=n
        neighbors = B(row(i)-1,col(i))+B(row(i),col(i)+1)+B(row(i)+1,col(i)+1)+B(row(i)-1,col(i)+1)+B(row(i)+1,col(i));
        if neighbors >=3 || neighbors ==1
            B(row(i),col(i))=0;
        else
            B(row(i),col(i))=1;
        end
    end

%Last column
elseif col(i)==n
    if row(i)~=1 && row(i)~=n
        neighbors = B(row(i)-1,col(i))+B(row(i),col(i)-1)+B(row(i)+1,col(i)-1)+B(row(i)-1,col(i)-1)+B(row(i)+1,col(i));
        if neighbors >=3 || neighbors ==1
            B(row(i),col(i))=0;
        else
            B(row(i),col(i))=1;
        end
    end
    
else 
    neighbors = B(row(i)-1,col(i)+1)+B(row(i)+1,col(i)+1)+B(row(i),col(i)+1)+B(row(i)-1,col(i))+B(row(i),col(i)-1)+B(row(i)+1,col(i)-1)+B(row(i)-1,col(i)-1)+B(row(i)+1,col(i));
        if neighbors >=3 || neighbors ==1
            B(row(i),col(i))=0;
        else
            B(row(i),col(i))=1;
        end
end
end
spy(B)
title(['Tempo=', num2str(t), ' de ', num2str(tf)])
pause(0.2)
end
figure
plot(pop);
end
