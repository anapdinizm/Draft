% function Draft_direction_search_min(p0,fun)
%função z=x.^2+y.^2, x=-4:0.1:4, y=-4:0.1:4

%Sejam sigma > 0, alpha e theta (0, 1) constantes dadas. Se x_k pertence IR^n é tal que Grad(f(x_k))~= 0, os passos para determinar x_{k+1} são:
%Passo 1: Escolher d_k pertencente IR^n, tal que
%(i) norm(d_k) >= norm(Grad f(x_k));
%(ii) Grad f'(x_k)d_k >= norm(?f(x_k)) norm(d_k).
%Passo 2: (Busca linear)
%(i) lambda = 1;
%(ii) Se f(x_k + lambda d_k) < f(x_k) + ???f'(x_k)d_k, ir a (iv);
%(iii) Escolher lambda_barra pertencente [0.1 lambda, 0.9lambda]. Fazer lambda = lambda_barra e ir a (ii);
%(iv) Fazer lambda_k = lambda, e x_{k+1} = x_k + lambda_k d_k.
v=-4:0.001:4;
[x,y]=meshgrid(v);
% if nargin<2
    fun=x.^2+y.^2;
    p0=[4,4];
% end
%constantes dadas
alpha=0.5;
theta=0.5;
sigma=0.3; %tamanho do passo, para que não vá rapidamente para zero

z=fun;
i=1;
p(i,:)=p0;
[rx,cx]=find(x==p0(1));
[ry,cy]=find(y==p0(2));
r=find(rx==ry);
c=find(cx==cy);
tol=1;
[gradx,grady] = gradient(z,.1,.1);
while (gradx(r,c)~=0 || grady(r,c)~=0) && tol>10^(-8)
    i
    %Passo 1
    norm_d=sigma*norm([gradx(r,c),grady(r,c)]);
    d(i,:)=[gradx(r,c),grady(r,c)]\(-theta*norm([gradx(r,c),grady(r,c)])*norm_d);
    %Passo 2
    lambda(i)=1;
    p_next=p(i,:)+lambda*d(i,:);
    %PRECISO VERIFICAR APROXIMAÇÃO
    [rx,cx]=find(x==p_next(1));
    [ry,cy]=find(y==p_next(2));
    r_next=find(rx==ry);
    c_next=find(cx==cy);
    %     if z(r_next,c_next)<z(r,c)+ alpha*lambda*[gradx(r,c),grady(r,c)]'*d(i,:)
    %         lambda(i+1)=lambda(i);
    %         p(i+1,:)=p_next;
    while z(r_next,c_next)>=z(r,c)+ alpha*lambda*[gradx(r,c),grady(r,c)]*d(i,:)
        lambda_barra=lambda-0.15;
        lambda=lambda_barra;
    end
    lambda(i+1)=lambda(i);
    p(i+1,:)=p_next;
    r=r_next;
    c=c_next;
    
    tol=p(i+1,:)-p(i,:)
    i=i+1;
    pause
end

   