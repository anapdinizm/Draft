%Simulação para o comportamento de alguns corpos que possuem atração
%gravitacional.
% n - # corpos

function Draft_Gravitacao(n)
% %Sistema de duas partículas
% m1= 3; %unidade da massa em Kg
% m2= 7;
% d = 1; %unidade de distância em m
% G = 6.673*10^-11; % unidade N*m^2/Kg
% P1=2; % posição do corpo 1
% P2=3; %posição do corpo 2
% v1=2;
% v2=3;
% delta_t=30; % tempo em segundos
% 
% F=G*m1*m2/d^2;
% % a1=F/m1;
% % a2=F/m2;
% 
% a1=(P2-P1)*G*m2/d^3;
% a2=(P1-P2)*G*m1/d^3;
% 
% %Supondo que o corpo 1 e 2 viajando em uma ...
%  ...tragetória retilínea com uma velocidade inicial v1 e v2
% v1_barra=v1+a1*delta_t;
% v2_barra=v2+a2*delta_t;
% 
% %Posição final do corpo em uma trajetória retilínea
% pos1=P1+v1_barra*delta_t;
% pos2=P2+v2_barra*delta_t;

%Força para n partículas
% unidade N*m^2/Kg
G = 6.673*10^-11; 
%m - unidade da massa em Kg
m=18*rand(n,1)+3;
Tl=tril(rand(n,n));
%d - unidade de distância em m
d=Tl*Tl';
d=d-diag(diag(d));
% v - unidade da velocidade metros/segundo
 v=12*rand(n,1);
%P_x, P_y - posições iniciais das particulas nas coordenadas...
... x e y
P_x=13*rand(n,1)+1;
P_y=15*rand(n,1)+0.5;
% tempo em segundos
delta_t=30; 

F=zeros(n,1);
for i=1:n
    for j=1:n
        if i~=j
            F(i)=F(i)+G*m(i)*m(j)/d(i,j)^2;
        end
    end
end
F

a=zeros(n,1);
for i=1:n
    a(i)=F(i)/m(i);
end
a

%Supondo que o corpo 1 e 2 estejam viajando em uma ...
 ...tragetória retilínea com uma velocidade inicial v1 e v2
for i=1:n
    v_barra(i)=v(i)+a(i)*delta_t;
end 
v_barra

%Posição final do corpo em uma trajetória retilínea
for i=1:n
    pos_x(i)=P_x(i)+v_barra(i)*delta_t;
    pos_y(i)=P_y(i)+v_barra(i)*delta_t;
end
pos_x
pos_y
figure(1)
plot(P_x,P_y,'.');
figure(2)
plot(pos_x,pos_y,'.');
end
