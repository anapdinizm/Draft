%Comportamento da equação de difusão em 1D
%Resolução pelo método implícito de Euler e por Crank-Nicolson.

function Draft_diffusion_1D
clear all
%% Entrada da dimensão em x
nx=50; %intervalos em x

%pontos em uma dimensão
nnx=nx;
nn=nnx+1;

L=2; %largura

dx=L/nx; %variação em x


%% Parâmetros para o tempo
nt=30; %número de intervalos em t
tf=20; %tempo final
dt=tf/nt; %variação em t

%% Parâmetros
D=0.125*10^-2; %coeficiente de viscosidade
u=0.003; %velocidade em x

%% Núcleo de péclet
disp('o núcleo de péclet é')
np=u*dx/D;
disp(np)
%% Expansão em Taylor
% Valores referentes as entradas da matriz M
% E= 1+(2*D *dt/(dx*dx)); %coef. de c_i
% F= (-(D/dx)+(u/2))*dt/dx; %coef. de c_{i+1}
% G= (-(D/dx)-(u/2))*dt/dx; %coef. de c_{i-1}

%% Crank-Nicolson matrizes
% Valores referentes as entradas da matriz M
E= 1+(D *dt/(dx*dx)); %coef. de c_i
F= (-(D/2*dx)+(u/4))*dt/dx; %coef. de c_{i+1}
G= (-(D/2*dx)-(u/4))*dt/dx; %coef. de c_{i-1}

% Valores referentes as entradas da matriz P
Ep= 1-(D *dt/(dx*dx)); %coef. de c_i
Fp= ((D/2*dx)-(u/4))*dt/dx; %coef. de c_{i+1}
Gp= ((D/2*dx)+(u/4))*dt/dx; %coef. de c_{i-1}

%% Construção da matriz M
M=E*eye(nn)+diag(diag(F*eye(nn-1)),1)+diag(diag(G*eye(nn-1)),-1);
%% Construção da matriz P
P=Ep*eye(nn)+diag(diag(Fp*eye(nn-1)),1)+diag(diag(Gp*eye(nn-1)),-1);

% Preciso incluir condição de Von neumann
%% Construção da imagem
C_ini=zeros(nn,1); %inicialização da matriz C_ini
C_ini(2,1)=2;

%dimensões em x e em y
dim_x=[0:dx:L];

% Rearranjando o vetor C na matriz conc
figure
for t=1:nt
    C=M\P*(C_ini); %Resolução do sistema
    C_ini=C;
    C_ini(nn,1)=0;  %Condição de von neumann
    C_ini(1,1)=0;
    
    h=plot(dim_x,C);
    xlabel('Dimensao de x')
    ylabel('Concentracao')
    title({['Difusão com D = ',num2str(D)];['u = ',num2str(u)];[ 'Núcleo de Péclet = ',num2str(np) ];['tempo (t) = ',num2str(t*dt),' de ',num2str(tf)]})
    shading interp
    
    pause(0.1)
end

end
