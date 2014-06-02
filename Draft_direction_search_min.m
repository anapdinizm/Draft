% Function try to find the minimum point of two variable function
%Algorithm Idea:
%Sejam sigma > 0, alpha e theta (0, 1) constantes dadas. Se x_k pertence IR^n é tal que Grad(f(x_k))~= 0, os passos para determinar x_{k+1} são:
%Passo 1: Escolher d_k pertencente IR^n, tal que
%(i) norm(d_k) >= norm(Grad f(x_k));
%(ii) Grad'(f(x_k))d_k >= norm(Grad(f(x_k))) norm(d_k).
%Passo 2: (Busca linear)
%(i) lambda = 1;
%(ii) Se f(x_k + lambda d_k) < f(x_k) + Grad'(f(x_k))d_k, ir a (iv);
%(iii) Escolher lambda_barra pertencente [0.1 lambda, 0.9lambda]. Fazer lambda = lambda_barra e ir a (ii);
%(iv) Fazer lambda_k = lambda, e x_{k+1} = x_k + lambda_k d_k.
function Draft_direction_search_min(p0,fun)

%Create the symbolic function
if nargin<2
    fun=x^3+y^2;
    p0=[1,0];
end

syms x y
% f(x,y)=x^3+y^2;
% f(x,y)=exp(x+2*y);
f(x,y)=fun;
gradx=diff(f,x,1);
grady=diff(f,y,1);
hessiana=[diff(gradx,x,1),diff(gradx,y,1);diff(grady,x,1),diff(grady,y,1)];


alpha=0.5;
theta=0.5;
sigma=0.3; %tamanho do passo, para que não vá rapidamente para zero

i=1;
p(i,:)=p0;
tol=[1,1];
lambda=1;

while (gradx(p0(1),p0(2))~=0 || grady(p0(1),p0(2))~=0) && norm(tol)>10^(-6)
    i
    %Passo 1
    norm_d=sigma*norm([gradx(p0(1),p0(2)),grady(p0(1),p0(2))]);
    norm_d=double(norm_d);
%     d0=[gradx(p0(1),p0(2)),grady(p0(1),p0(2))]\(-theta*norm([gradx(p0(1),p0(2)),grady(p0(1),p0(2))])*norm_d);
    d0=-theta*[gradx(p0(1),p0(2)),grady(p0(1),p0(2))];
    d0=double(d0);
    %Passo 2
    p_next=p0+lambda*d0;
    lambda_barra=lambda;
    while (lambda_barra > 0.1*lambda) && (f(p_next(1),p_next(2)) >= f(p0(1),p0(2))+ alpha*lambda_barra*[gradx(p0(1),p0(2)),grady(p0(1),p0(2))]*d0')
        lambda_barra=lambda-0.15;
    end
    lambda=lambda_barra;
    p(i+1,:)=p_next;
    tol=p_next-p0;
    p0=p_next ;
    d(i,:)=d0;
    i=i+1;
    pause
end
   
% %plot
v=-4:0.1:4;
[x,y]=meshgrid(v);
fun=f(x,y);
fun=double(fun);
[Ngradx,Ngrady] = gradient(fun,.1,.1);
contour(x,y,fun);
hold on
quiver(x,y,Ngradx,Ngrady);
scatter(p(:,1),p(:,2),100,'g','fill')
plot(p(:,1),p(:,2))
hold off
figure
surf(x,y,fun)

