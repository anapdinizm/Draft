%Comportamento da equa��o de difus�o em 2D
%Resolu��o pelo m�todo impl�cito de Euler .

function Draft_diffusion_2D_euler_backward
clear all
%% Entradas das dimens�es da malha
nx=25; %intervalos em x
ny=5; %intervalos em y

%pontos na malha
nnx=nx;
nny=ny+1;
nn=nnx*nny;

L=2; %largura
H=0.5; %altura

dx=L/nx; %varia��o em x
dy=H/ny; %varia��o em y

%% Par�metros para o tempo
nt=200; %n�mero de intervalos em t
tf=20; %tempo final
dt=tf/nt; %varia��o em t

%% Par�metros
D=0.125*10^-2; %coeficiente de viscosidade
u=0.01; %velocidade em x
v=0.005; %velocidade em y
mu=0.05; %fator de decaimento

%% Defini��o do vetor Fonte
f=zeros(nn,1);
%  f(nny,1)=0.2412;
%  f(20*nny+2,1)=0.2412;
% f(20*nny+6,1)=0.2412;
% f(20*nny+5,1)=10;
% f(19*nny+5,1)=0.2412;
% f(19*nny+6,1)=0.2412;
% f(21*nny+5,1)=0.2412;
% f(21*nny+6,1)=0.2412;

%% N�cleo de p�clet
disp('o n�cleo de p�clet �')
np=v*dy/D;
disp(np)

%% Euler Impl�cito
%Valores referentes as entradas da matriz M
A= (-(D/dx)+(u/2))*dt/dx; %coef. de c_{i+nny}
B= -((D/dx)+(u/2))*dt/dx; %coef. de c_{i-nny}
E= 1+((2*D/(dx*dx))+(2*D/(dy*dy))+ mu)*dt; %coef. de c_i
F= (-(D/dy)+(v/2))*dt/dy; %coef. de c_{i+1}
G= -((D/dy)+(v/2))*dt/dy; %coef. de c_{i-1}


%% Constru��o da matriz M
M=E*eye(nn)+diag(diag(F*eye(nn-1)),1)+diag(diag(G*eye(nn-1)),-1)+diag(diag(B*eye(nn-nny)),-nny)+diag(diag(A*eye(nn-nny)),nny);


%% Condi��o das bordas inferior, superior e lateral direita 
for i=1:nn
    if mod(i,nny)==1 %condicao para a borda inferior
        if i==1 % canto esquerdo
            M(i,i+1)=0;%F+G;
        else %demais condi��es para a borda inferior
            M(i,i-1)=0;
            M(i,i-1)=0;
            M(i,i+1)=F+G; %condicao -2*D*nt/ny^2
            if i==nn-nny+1 % canto direito
                M(i,i-nny)=0;%A+B; %condicao -2*D*nt/nx^2
            end
        end
        
        
    elseif mod(i,nny)==0 % condicao para a borda superior
        if i==nny % canto direito
            M(i,i+1)=0;
            M(i,i-1)=0;%F+G; %condicao -2*D*nt/ny^2
        else %demais condi��es para a borda superior
            if i==nn % canto direito
                M(i,i-nny)=0;%A+B; %condicao -2*D*nt/nx^2
            else
                M(i,i+1)=0;
                M(i,i-1)=F+G; %condicao -2*D*nt/ny^2
            end
        end
        
    elseif i>nn-nny+1 && i<nn %condi��o da borda direita
        M(i,i-nny)=0;%A+B; %condicao -2*D*nt/nx^2
        
    end
end

%% Constru��o da imagem
C_ini=zeros(nn,1); %inicializa��o da matriz C_ini
C_ini((nny+4),1)=10; 
% conc=zeros(nnx+1,nny); %inicializa��o da matriz b

%dimens�es em x e em y
dim_x=[0:dx:L];
dim_y=[0:dy:H];

% Rearranjando o vetor C na matriz conc
figure
for t=1:nt
    C=M\(C_ini+dt*f); %Resolu��o do sistema
    C_ini=C;
%     for k=(nn-nny+1):nn
%         C_ini(k,1)=0; %condi��o de Von Neumann
%     end
%     for k=1:nny-1
%         C_ini(k,1)=0;
%     end
    
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
    title({['Difus�o com D = ',num2str(D)];['u = ',num2str(u) ', v = ',num2str(v) ', {\mu} = ',num2str(mu)];[ 'N�cleo de P�clet= ',num2str(vn) ];['tempo (t) = ',num2str(t*dt),' de ',num2str(tf)]})
    xlabel('Coordenada em (x)')
    ylabel('Coordenada em (y)')
    zlabel('Perfil da propriedade de transporte')
    if t==1
        saveas(h, 'fig1_p1_e6_2d', 'fig');
        saveas(h, 'fig1_p1_e6_2d', 'jpg');
        saveas(h, 'fig1_p1_e6_2d', 'pdf');
        pause;
    end
%     if t==(nt/4)
%         saveas(h, 'fig3_fonte5', 'fig');
%         saveas(h, 'fig3_fonte5', 'jpg');
%         saveas(h, 'fig3_fonte5', 'pdf');
%         pause;
%     end
%     if t==(nt/2)
%         saveas(h, 'fig5_fonte5', 'fig');
%         saveas(h, 'fig5_fonte5', 'jpg');
%         saveas(h, 'fig5_fonte5', 'pdf');
%         pause;
%     end
    drawnow;
    refreshdata(h)
saveas(h, 'fig4_p1_e6_2d', 'fig');
saveas(h, 'fig4_p1_e6_2d', 'jpg');
saveas(h, 'fig4_p1_e6_2d', 'pdf');
end
end
