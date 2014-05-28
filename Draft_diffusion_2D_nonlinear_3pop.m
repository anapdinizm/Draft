%Funcao nao linear de dispersao da interacao entre tres populacoes com crescimento de Verhulst
%Soluco para a dependencia da matriz pentadiagonal do valor inicial da
%propria população com o metodo "preditor corretor"
%Exemplo cobra - gavi�o - roedor
function Draft_diffusion_2D_nonlinear_3pop
clear all
%% Entradas das dimensoes da malha
nx=40; %intervalos em x
ny=10; %intervalos em y

%pontos na malha
nnx=nx;
nny=ny+1;
nn=nnx*nny;

L=2; %largura
H=0.5; %altura

dx=L/nx; %variacao em x
dy=H/ny; %variacao em y

%% Parametros para o tempo
nt=200; %numero de intervalos em t
tf=20; %tempo final
dt=tf/nt; %variacao em t

%% Parametros
D=0.125*10^-2; %coeficiente de viscosidade
u=0; %velocidade em x
v=0.001; %velocidade em y
lambda=0.2;
gama=0.05;
mu=gama-lambda; %fator de decaimento

K=20; %capacidade de suporte

%% Definição do vetor Fonte
f=zeros(nn,1);
% f(20*nny+6,1)=0.2412;

%% N�cleo de p�clet
disp('o núcleo de péclet é')
np=v*dy/(D);
disp(np)

%% Euler Impl�cito
%Valores referentes as entradas da matriz M
A= (-(D/dx)+(u/2))*dt/dx; %coef. de c_{i+nny}
B= -((D/dx)+(u/2))*dt/dx; %coef. de c_{i-nny}
F= (-(D/dy)+(v/2))*dt/dy; %coef. de c_{i+1}
G= -((D/dy)+(v/2))*dt/dy; %coef. de c_{i-1}

%%Condição INICIAL***
C0=zeros(nn,1);
C0(55,1)=1;
Caux=C0;
Caux1=C0;

conc=zeros(nnx+1,nny);

dim_x=[0:dx:L];
dim_y=[0:dy:H];

figure
for t=1:nt
    
    %%Montando a matriz do sistema com a condição inicial
    %Encontrando a melhor aproximação para o valor inicial, "método preditor corretor" loop das estrelas
    for estrela=1:7
        %coeficiente relativo à diagonal principal
        E= 1+((2*D/(dx*dx))+(2*D/(dy*dy))+ mu +(lambda/K)*Caux)*dt ;
        M=diag(E)+diag(G*ones(nn-1,1),-1)+diag(F*ones(nn-1,1),1)+diag(B*ones(nn-nny,1),-nny)+diag(A*ones(nn-nny,1),nny);
        
        %Arrumando elementos referentes às condições de contorno
        
        %fronteira de baixo
        for k=nny+1:nny:nn-nny
            M(k,k-1)=0;
            M(k,k+1)=-2*D*dt/(dy*dy);
        end
        %fronteira de cima
        for k=nny:nny:nn-nny
            M(k,k+1)=0;
            M(k,k-1)=-2*D*dt/(dy*dy);
        end
        %fronteira de direita
        for k=(nnx-1)*nny+2:nn-1
            M(k,k-nny)=-2*D*dt/(dx*dx);
        end
        
        %relativo à primeira equação
        M(1,2)=-2*D*dt/(dy*dy);
        
        %relativo a equação do canto superior esquerdo
        M(nny,nny-1)=-2*D*dt/(dy*dy);
        M(nny,nny+1)=0;
        
        %relativo a equação do canto inferior direito
        M(nn-nny+1,nn-nny)=0;
        M(nn-nny+1,nn-2*nny+1)=-2*D*dt/(dx*dx);
        M(nn-nny+1,nn-nny+2)=-2*D*dt/(dy*dy);
        
        %relativo a equação do canto superior direito
        M(nn, nn-1)=-2*D*dt/(dy*dy);
        M(nn, nn-nny)=-2*D*dt/(dx*dx);
        
        
        Caux1=M\(C0+dt*f);
        Caux=Caux1;
        %***
    end
    
    %****
    
    E= 1+((2*D/(dx*dx))+(2*D/(dy*dy))+ mu +(lambda/K)*Caux)*dt ;
    M=diag(E)+diag(G*ones(nn-1,1),-1)+diag(F*ones(nn-1,1),1)+diag(B*ones(nn-nny,1),-nny)+diag(A*ones(nn-nny,1),nny);
    
    %Arrumando elementos referentes às condições de contorno
    
    %fronteira de baixo
    for k=nny+1:nny:nn-nny
        M(k,k-1)=0;
        M(k,k+1)=-2*D*dt/(dy*dy);
    end
    %fronteira de cima
    for k=nny:nny:nn-nny
        M(k,k+1)=0;
        M(k,k-1)=-2*D*dt/(dy*dy);
    end
    %fronteira de direita
    for k=(nnx-1)*nny+2:nn-1
        M(k,k-nny)=-2*D*dt/(dx*dx);
    end
    
    %relativo à primeira equação
    M(1,2)=-2*D*dt/(dy*dy);
    
    %relativo a equação do canto superior esquerdo
    M(nny,nny-1)=-2*D*dt/(dy*dy);
    M(nny,nny+1)=0;
    
    %relativo a equação do canto inferior direito
    M(nn-nny+1,nn-nny)=0;
    M(nn-nny+1,nn-2*nny+1)=-2*D*dt/(dx*dx);
    M(nn-nny+1,nn-nny+2)=-2*D*dt/(dy*dy);
    
    %relativo a equação do canto superior direito
    M(nn, nn-1)=-2*D*dt/(dy*dy);
    M(nn, nn-nny)=-2*D*dt/(dx*dx);
    
    
    C=M\(C0+dt*f);
    %****
    C=M\(C0+dt*f);
    C0=C;
    conc=zeros(nny,nnx+1);
    for k=2:nnx+1
        I=(k-1)*nny;
        for i=1:nny
            conc(i,k)=C(I,1);
            I=I-1;
        end
    end
    %Plot da superf�cie
    h=surf(dim_x,dim_y,conc,'EdgeColor','none');
    %     shading interp
    colorbar();
    title({['Difusão com D = ',num2str(D)];['u = ',num2str(u) ', v = ',num2str(v) ', {\mu} = ',num2str(mu)];[ 'Núcleo de Péclet= ',num2str(np) ];['tempo (t) = ',num2str(t*dt),' de ',num2str(tf)]})
    xlabel('Coordenada em (x)')
    ylabel('Coordenada em (y)')
    zlabel('Perfil da propriedade de transporte')
    drawnow;
    refreshdata(h)
end
end
