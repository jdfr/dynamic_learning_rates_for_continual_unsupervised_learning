function [Samples]=GenerateManifolds3D(NdxDataset,Precision,PlotIt)
% Generar variedades (manifolds) de algunos conjuntos de datos 3D
% NdxDataset=Número del conjunto de datos que se desea (1 a 9)
% Precision=Precisión del manifold
% PlotIt=Flag to generate a plot of the manifold
% Tomado de Todd Wittman, MANIfold Learning Matlab Demo, http://www.math.ucla.edu/~wittman/mani/index.html
% Añadida la Ueda's spiral


% ExParam=Parámetro extra (altura)
N=Precision;
ExParam=1;
switch NdxDataset
    case 1  % Swiss Roll
        % Muestras=GenerarManifolds3D(1,20);
  
        tt = (3*pi/2)*(1+2*(0:(1/(2*N)):1));  
        height = 21*(0:(1/N):1);
        

        xx = (tt.*cos(tt))'*ones(size(height));
        yy = ones(size(tt))'*height;
        zz = (tt.*sin(tt))'*ones(size(height));
        cc = tt'*ones(size(height));
        if PlotIt
            for Ndx1=1:(size(xx,1)-1)
                for Ndx2=1:(size(xx,2)-1)
                    patch([xx(Ndx1,Ndx2) xx(Ndx1+1,Ndx2) xx(Ndx1+1,Ndx2+1) xx(Ndx1,Ndx2+1) ],...
                            [yy(Ndx1,Ndx2) yy(Ndx1+1,Ndx2) yy(Ndx1+1,Ndx2+1) yy(Ndx1,Ndx2+1) ],...
                            [zz(Ndx1,Ndx2) zz(Ndx1+1,Ndx2) zz(Ndx1+1,Ndx2+1) zz(Ndx1,Ndx2+1) ],...
                            [cc(Ndx1,Ndx2) cc(Ndx1+1,Ndx2) cc(Ndx1+1,Ndx2+1) cc(Ndx1,Ndx2+1) ]);
                end
            end
        end
        Samples=[xx(:)';yy(:)';zz(:)'];

    case 2  % Swiss Hole
        %Muestras=GenerarManifolds3D(2,20);
        % Swiss Roll w/ hole example taken from Donoho & Grimes
        tt = (3*pi/2)*(1+2*(0:(1/(2*N)):1));  
        height = 21*(0:(1/N):1);
        kl=zeros(numel(tt),numel(height));
        for ii = 1:numel(tt)
            for jj=1:numel(height)
                if ( (tt(ii) > 9) && (tt(ii) < 12))
                    if ((height(jj) > 9) && (height(jj) <14))
                        kl(ii,jj) = 1;
                    end;
                end;
            end
        end;
        
        xx = (tt.*cos(tt))'*ones(size(height));
        yy = ones(size(tt))'*height;
        zz = (tt.*sin(tt))'*ones(size(height));
        cc = tt'*ones(size(height));
        if PlotIt
            for Ndx1=1:(size(xx,1)-1)
                for Ndx2=1:(size(xx,2)-1)
                    if (kl(Ndx1,Ndx2)==0) && (kl(Ndx1+1,Ndx2)==0) && (kl(Ndx1,Ndx2+1)==0) && ...
                            (kl(Ndx1+1,Ndx2+1)==0)
                        patch([xx(Ndx1,Ndx2) xx(Ndx1+1,Ndx2) xx(Ndx1+1,Ndx2+1) xx(Ndx1,Ndx2+1) ],...
                            [yy(Ndx1,Ndx2) yy(Ndx1+1,Ndx2) yy(Ndx1+1,Ndx2+1) yy(Ndx1,Ndx2+1) ],...
                            [zz(Ndx1,Ndx2) zz(Ndx1+1,Ndx2) zz(Ndx1+1,Ndx2+1) zz(Ndx1,Ndx2+1) ],...
                            [cc(Ndx1,Ndx2) cc(Ndx1+1,Ndx2) cc(Ndx1+1,Ndx2+1) cc(Ndx1,Ndx2+1) ]);
                    end
                end
            end
        end
        Samples=[xx(kl(:)'==0);yy(kl(:)'==0);zz(kl(:)'==0)];
        
    case 3  % Corner Planes
        % Muestras=GenerarManifolds3D(3,100);
        k = 1;
        MyParam=1000;
        xMax = floor(sqrt(MyParam));
        yMax = ceil(MyParam/xMax);
        cornerPoint = floor(yMax/2);
        for x = 0:xMax
            for y = 0:yMax
                if y <= cornerPoint
                    xx(x+1,y+1)=x;
                    yy(x+1,y+1)=y;
                    zz(x+1,y+1)=0;
                    cc(x+1,y+1)=y;
                    X(k,:) = [x,y,0];
                    ColorVector(k) = y;
                else
                    xx(x+1,y+1)=x;
                    yy(x+1,y+1)=cornerPoint+(y-cornerPoint)*cos(pi*ExParam/180);
                    zz(x+1,y+1)=(y-cornerPoint)*sin(pi*ExParam/180);     
                    cc(x+1,y+1)=y;
                    X(k,:) = [x,cornerPoint+(y-cornerPoint)*cos(pi*ExParam/180),(y-cornerPoint)*sin(pi*ExParam/180)];
                    ColorVector(k) = y;
                end;
                k = k+1;
            end;
        end;

        if PlotIt
            for Ndx1=1:(size(xx,1)-1)
                for Ndx2=1:(size(xx,2)-1)
                    patch([xx(Ndx1,Ndx2) xx(Ndx1+1,Ndx2) xx(Ndx1+1,Ndx2+1) xx(Ndx1,Ndx2+1) ],...
                            [yy(Ndx1,Ndx2) yy(Ndx1+1,Ndx2) yy(Ndx1+1,Ndx2+1) yy(Ndx1,Ndx2+1) ],...
                            [zz(Ndx1,Ndx2) zz(Ndx1+1,Ndx2) zz(Ndx1+1,Ndx2+1) zz(Ndx1,Ndx2+1) ],...
                            [cc(Ndx1,Ndx2) cc(Ndx1+1,Ndx2) cc(Ndx1+1,Ndx2+1) cc(Ndx1,Ndx2+1) ]);
                end
            end
        end

        Samples=[xx(:)';yy(:)';zz(:)'];

    case 4  % Punctured Sphere by Saul & Roweis
        % Muestras=GenerarManifolds3D(4,1000);
        inc = 9/sqrt(N); 
        [xx,yy] = meshgrid(-5:inc:5);
        rr2 = xx(:).^2 + yy(:).^2;
        [tmp ii] = sort(rr2);
        a = 4./(4+(xx.^2+yy.^2));
        xx=a.*xx;
        yy=a.*yy;
        zz=ExParam*2*(1-a);
        cc=zz;
        if PlotIt
            for Ndx1=1:(size(xx,1)-1)
                for Ndx2=1:(size(xx,2)-1)
                    patch([xx(Ndx1,Ndx2) xx(Ndx1+1,Ndx2) xx(Ndx1+1,Ndx2+1) xx(Ndx1,Ndx2+1) ],...
                        [yy(Ndx1,Ndx2) yy(Ndx1+1,Ndx2) yy(Ndx1+1,Ndx2+1) yy(Ndx1,Ndx2+1) ],...
                        [zz(Ndx1,Ndx2) zz(Ndx1+1,Ndx2) zz(Ndx1+1,Ndx2+1) zz(Ndx1,Ndx2+1) ],...
                        [cc(Ndx1,Ndx2) cc(Ndx1+1,Ndx2) cc(Ndx1+1,Ndx2+1) cc(Ndx1,Ndx2+1) ]);
                end
            end
        end
        Samples=[xx(:)';yy(:)';zz(:)'];

    case 5  % Twin Peaks by Saul & Roweis
        %Muestras=GenerarManifolds3D(5,200);
        inc = 1.5 / sqrt(N); 
        [xx,yy] = meshgrid(-1:inc:1);
        zz = sin(pi*xx).*tanh(3*yy);
        cc=zz;
        if PlotIt
            for Ndx1=1:(size(xx,1)-1)
                for Ndx2=1:(size(xx,2)-1)
                    patch([xx(Ndx1,Ndx2) xx(Ndx1+1,Ndx2) xx(Ndx1+1,Ndx2+1) xx(Ndx1,Ndx2+1) ],...
                            [yy(Ndx1,Ndx2) yy(Ndx1+1,Ndx2) yy(Ndx1+1,Ndx2+1) yy(Ndx1,Ndx2+1) ],...
                            [zz(Ndx1,Ndx2) zz(Ndx1+1,Ndx2) zz(Ndx1+1,Ndx2+1) zz(Ndx1,Ndx2+1) ],...
                            [cc(Ndx1,Ndx2) cc(Ndx1+1,Ndx2) cc(Ndx1+1,Ndx2+1) cc(Ndx1,Ndx2+1) ]);
                end
            end
        end
        Samples=[xx(:)';yy(:)';zz(:)'];

    case 6  % Toroidal Helix by Coifman & Lafon
        % Muestras=GenerarManifolds3D(7,500);
        t = (1:N)'/N;
        t = t.^(ExParam)*2*pi;
        Samples = [(2+cos(8*t)).*cos(t) (2+cos(8*t)).*sin(t) sin(8*t)];
        if PlotIt
            MapaColor=colormap;
            NumColores=size(MapaColor,1);
            for Ndx=1:(size(Samples,1)-1)
                MiHandle=line([Samples(Ndx,1) Samples(Ndx+1,1)],[Samples(Ndx,2) Samples(Ndx+1,2)],...
                    [Samples(Ndx,3) Samples(Ndx+1,3)]);
                set(MiHandle,'Color',MapaColor(ceil((1+Samples(Ndx,3))*0.5*NumColores),:))
            end
        end
        Samples=Samples';

    case 7  % Ueda's spiral
        % Muestras=GenerarManifolds3D(9,500);
        t = (((0:N-1)/(N-1))*4*pi)';
        Samples = [(13-(0.5*t)).*cos(t) -(13-(0.5*t)).*sin(t) t] ;        
        MapaColor=colormap;
        NumColores=size(MapaColor,1);
        Maximo=max(t);
        for Ndx=1:(size(Samples,1)-1)
            MiHandle=line([Samples(Ndx,1) Samples(Ndx+1,1)],[Samples(Ndx,2) Samples(Ndx+1,2)],...
                [Samples(Ndx,3) Samples(Ndx+1,3)]);
            set(MiHandle,'Color',MapaColor(1+ceil((Samples(Ndx,3)/Maximo)*(NumColores-1)),:))
        end
        Samples=Samples';
end

