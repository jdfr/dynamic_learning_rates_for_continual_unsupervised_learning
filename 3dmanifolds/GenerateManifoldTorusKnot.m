function [Muestras]=GenerateManifoldTorusKnot(p,q,Precision)
% Precision=Cantidad de facetas. Por ejemplo, 100.


Phi=2*p*pi*(0:(1/Precision):1);

% Puntos de la curva
Centros(1,:)=(2+cos(q*Phi/p)).*cos(Phi);
Centros(2,:)=(2+cos(q*Phi/p)).*sin(Phi);
Centros(3,:)=sin(q*Phi/p);

% Vectores directores en esos puntos de la curva, (dx / d phi, dy / d Phi, dz
% /d Phi ). Obtenidos con sym/diff
Direc(1,:)=-sin(Phi).*(cos((Phi*q)/p) + 2) - (q*sin((Phi*q)/p).*cos(Phi))/p;
Direc(2,:)=cos(Phi).*(cos((Phi*q)/p) + 2) - (q*sin((Phi*q)/p).*sin(Phi))/p;
Direc(3,:)=(q*cos((Phi*q)/p))/p;

% Norma del vector director, hallada con sym/simple
%MyNorm=(cos((Phi*q)/p)^2 + 8*cos((Phi*q)/(2*p))^2 + q^2/p^2)^(1/2);


% Angulos a considerar 
NumAngulos=8;
Theta=0:(2*pi/NumAngulos):(2*pi);

% Vectores normales a la curva en esos puntos
NumColores=size(colormap,1);
MapaColores=colormap;
for Ndx=1:size(Direc,2)
    MiNormA(:,Ndx)=cross(Direc(:,Ndx),[1 1 1]);
    MiNormB(:,Ndx)=cross(Direc(:,Ndx),MiNormA(:,Ndx)); 
    
    if Ndx<0.5*size(Direc,2)
        NdxColor=ceil(2*Ndx*NumColores/size(Direc,2));
    else
        NdxColor=NumColores-ceil((Ndx-0.5*size(Direc,2))*2*NumColores/size(Direc,2));
    end    
        
    for Ndx2=1:numel(Theta)
        MiTheta=Theta(Ndx2);
        Delta=cos(MiTheta)*MiNormA(:,Ndx)+sin(MiTheta)*MiNormB(:,Ndx);
        Delta=0.1*Delta/norm(Delta);
        xx(Ndx,Ndx2)=Centros(1,Ndx)+Delta(1);
        yy(Ndx,Ndx2)=Centros(2,Ndx)+Delta(2);
        zz(Ndx,Ndx2)=Centros(3,Ndx)+Delta(3);
        cc(Ndx,Ndx2)=NdxColor;
    end
end

for Ndx1=1:(size(xx,1)-1)
    for Ndx2=1:(size(xx,2)-1)
        patch([xx(Ndx1,Ndx2) xx(Ndx1+1,Ndx2) xx(Ndx1+1,Ndx2+1) xx(Ndx1,Ndx2+1) ],...
                [yy(Ndx1,Ndx2) yy(Ndx1+1,Ndx2) yy(Ndx1+1,Ndx2+1) yy(Ndx1,Ndx2+1) ],...
                [zz(Ndx1,Ndx2) zz(Ndx1+1,Ndx2) zz(Ndx1+1,Ndx2+1) zz(Ndx1,Ndx2+1) ],...
                [cc(Ndx1,Ndx2) cc(Ndx1+1,Ndx2) cc(Ndx1+1,Ndx2+1) cc(Ndx1,Ndx2+1) ]);
            %                'EdgeColor',MapaColores(cc(Ndx1,Ndx2),:));
    end
end
Muestras=[xx(:);yy(:);zz(:)];
