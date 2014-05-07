%Comportamento da equa��o de difus�o em 1D
%Resolu��o pelo m�todo impl�cito de Euler e por Crank-Nicolson.

function Draft_diffusion_1D
clear all
%% Entrada da dimens�o em x
nx=50; %intervalos em x

%pontos em uma dimens�o
nnx=nx;
nn=nnx+1;

L=2; %largura

dx=L/nx; %varia��o em x


%% Par�metros para o tempo
nt=30; %n�mero de intervalos em t
tf=20; %tempo final
dt=tf/nt; %varia��o em t

%% Par�metros
D=0.125*10^-2; %coeficiente de viscosidade
u=0.003; %velocidade em x

%% N�cleo de p�clet
disp('o n�cleo de p�clet �')
np=u*dx/D;
disp(np)
%% Expans�o em Taylor
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

%% Constru��o da matriz M
M=E*eye(nn)+diag(diag(F*eye(nn-1)),1)+diag(diag(G*eye(nn-1)),-1);
%% Constru��o da matriz P
P=Ep*eye(nn)+diag(diag(Fp*eye(nn-1)),1)+diag(diag(Gp*eye(nn-1)),-1);

% Preciso incluir condi��o de Von neumann
%% Constru��o da imagem
C_ini=zeros(nn,1); %inicializa��o da matriz C_ini
C_ini(2,1)=2;

%dimens�es em x e em y
dim_x=[0:dx:L];

% Rearranjando o vetor C na matriz conc
figure
for t=1:nt
    C=M\P*(C_ini); %Resolu��o do sistema
    C_ini=C;
    C_ini(nn,1)=0;  %Condi��o de von neumann
    C_ini(1,1)=0;
    
    h=plot(dim_x,C);
    xlabel('Dimensao de x')
    ylabel('Concentracao')
    title({['Difus�o com D = ',num2str(D)];['u = ',num2str(u)];[ 'N�cleo de P�clet = ',num2str(np) ];['tempo (t) = ',num2str(t*dt),' de ',num2str(tf)]})
    shading interp
    
    pause(0.1)
end

end
