%Comportamento da equação de difusão em 1D
%Resolução pelo método implícito de Euler e por Crank-Nicolson.

function Draft_diffusion_2D
clear all
%% Entradas das dimensões da malha
nx=40; %intervalos em x
ny=10; %intervalos em y

%pontos na malha
nnx=nx;
nny=ny+1;
nn=nnx*nny;

L=2; %largura
H=0.5; %altura

dx=L/nx; %variação em x
dy=H/ny; %variação em y

%% Parâmetros para o tempo
nt=200; %número de intervalos em t
tf=20; %tempo final
dt=tf/nt; %variação em t

%% Parâmetros
D=0.125*10^-2; %coeficiente de viscosidade
u=0.01; %velocidade em x
v=0.005; %velocidade em y
mu=0.05; %fator de decaimento

%% Definição do vetor Fonte
f=zeros(nn,1);
%  f(nny,1)=0.2412;
%  f(20*nny+2,1)=0.2412;
% f(20*nny+6,1)=0.2412;
% f(20*nny+5,1)=10;
% f(19*nny+5,1)=0.2412;
% f(19*nny+6,1)=0.2412;
% f(21*nny+5,1)=0.2412;
% f(21*nny+6,1)=0.2412;

%% Núcleo de péclet
disp('o núcleo de péclet é')
np=v*dy/D;
disp(np)

%% Expansão em Taylor
%Valores referentes as entradas da matriz M
% A= (-(D/dx)+(u/2))*dt/dx; %coef. de c_{i+nny}
% B= -((D/dx)+(u/2))*dt/dx; %coef. de c_{i-nny}
% E= 1+((2*D/(dx*dx))+(2*D/(dy*dy))+ mu)*dt; %coef. de c_i
% F= (-(D/dy)+(v/2))*dt/dy; %coef. de c_{i+1}
% G= -((D/dy)+(v/2))*dt/dy; %coef. de c_{i-1}

%% Crank-Nicolson matrizes
% Valores referentes as entradas da matriz M
A= (-(D/2*dx)+(u/4))*dt/dx; %coef. de c_{i+nny}
B= (-(D/2*dx)-(u/4))*dt/dx; %coef. de c_{i-nny}
E= 1+(D *dt/(dx*dx))+(D *dt/(dy*dy)); %coef. de c_i
F= (-(D/2*dy)+(v/4))*dt/dy; %coef. de c_{i+1}
G= (-(D/2*dy)-(v/4))*dt/dy; %coef. de c_{i-1}

% Valores referentes as entradas da matriz P
Ap= ((D/2*dx)-(u/4))*dt/dx; %coef. de c_{i+nny}
Bp= ((D/2*dx)+(u/4))*dt/dx; %coef. de c_{i-nny}
Ep= 1-(D *dt/(dx*dx))-(D *dt/(dy*dy)); %coef. de c_i
Fp= ((D/2*dy)-(v/4))*dt/dy; %coef. de c_{i+1}
Gp= ((D/2*dy)+(v/4))*dt/dy; %coef. de c_{i-1}



%% Construção da matriz M
M=E*eye(nn)+diag(diag(F*eye(nn-1)),1)+diag(diag(G*eye(nn-1)),-1)+diag(diag(B*eye(nn-nny)),-nny)+diag(diag(A*eye(nn-nny)),nny);

%% Construção da matriz P
P=Ep*eye(nn)+diag(diag(Fp*eye(nn-1)),1)+diag(diag(Gp*eye(nn-1)),-1);

%% Condição das bordas inferior, superior e lateral direita 
for i=1:nn
    if mod(i,nny)==1 %condicao para a borda inferior
        if i==1 % canto esquerdo
            M(i,i+1)=F+G;
        else %demais condições para a borda inferior
            M(i,i-1)=0;
            M(i,i+1)=F+G; %condicao -2*D*nt/ny^2
            if i==nn-nny+1 % canto direito
                M(i,i-nny)=A+B; %condicao -2*D*nt/nx^2
            end
        end
        
        
    elseif mod(i,nny)==0 % condicao para a borda superior
        if i==nny % canto esquerdo
            M(i,i+1)=0;
            M(i,i-1)=F+G; %condicao -2*D*nt/ny^2
        else %demais condições para a borda superior
            if i==nn % canto direito
                M(i,i-nny)=A+B; %condicao -2*D*nt/nx^2
            else
                M(i,i+1)=0;
                M(i,i-1)=F+G; %condicao -2*D*nt/ny^2
            end
        end
        
    elseif i>nn-nny+1 && i<nn %condição da borda esquerda
        M(i,i-nny)=A+B; %condicao -2*D*nt/nx^2
        
    end
end

%% Construção da imagem
C_ini=zeros(nn,1); %inicialização da matriz C_ini
C_ini(nny+6,1)=10; 
% conc=zeros(nnx+1,nny); %inicialização da matriz b

%dimensões em x e em y
dim_x=[0:dx:L];
dim_y=[0:dy:H];

% Rearranjando o vetor C na matriz conc
figure
for t=1:nt
    C=M\P*(C_ini+dt*f); %Resolução do sistema
    C_ini=C;
    C_ini(nn-nny+1:nn,1)=0; %condição de Von Neumann
    C_ini(1:nny-1,1)=0;
    
    conc=zeros(nny,nnx+1);
    for k=2:nnx+1
        I=(k-1)*nny;
        for i=1:nny
            conc(i,k)=C(I,1);
            I=I-1;
        end
    end
    
    %Plot da superfície
    h=surf(dim_x,dim_y,conc,'EdgeColor','none');
%     shading interp
    colorbar();
    title({['Difusão com D = ',num2str(D)];['u = ',num2str(u) ', v = ',num2str(v) ', {\mu} = ',num2str(mu)];['tempo (t) = ',num2str(t*dt),' de ',num2str(tf)]})
    xlabel('Coordenada em (x)')
    ylabel('Coordenada em (y)')
    zlabel('Perfil da propriedade de transporte')
    if t==1
%         saveas(h, 'fig1_fonte5', 'fig');
%         saveas(h, 'fig1_fonte5', 'jpg');
%         saveas(h, 'fig1_fonte5', 'pdf');
        pause;
    end
    if t==(nt/4)
%         saveas(h, 'fig3_fonte5', 'fig');
%         saveas(h, 'fig3_fonte5', 'jpg');
%         saveas(h, 'fig3_fonte5', 'pdf');
        pause;
    end
    if t==(nt/2)
%         saveas(h, 'fig5_fonte5', 'fig');
%         saveas(h, 'fig5_fonte5', 'jpg');
%         saveas(h, 'fig5_fonte5', 'pdf');
        pause;
    end
    drawnow;
    refreshdata(h)
% saveas(h, 'fig7_fonte5', 'fig');
% saveas(h, 'fig7_fonte5', 'jpg');
% saveas(h, 'fig7_fonte5', 'pdf');
end
end