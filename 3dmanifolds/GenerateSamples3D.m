function [Samples]=GenerateSamples3D(NdxDataset,NumSamples,StdNoise)
% Generar muestras de algunos conjuntos de datos 3D
% NdxDataset=Número del conjunto de datos que se desea (1 a 9)
% NumMuestras=Número de muestras que generar
% StdNoise=Standard deviation of the spherical Gaussian noise
% Tomado de Todd Wittman, MANIfold Learning Matlab Demo, http://www.math.ucla.edu/~wittman/mani/index.html
% Añadida la Ueda's spiral


% ExParam=Parámetro extra (altura)
N=NumSamples;
ExParam=1;
switch NdxDataset
    case 1  % Swiss Roll
        tt = (3*pi/2)*(1+2*rand(1,N));  
        height = 21*rand(1,N);
        Samples = [tt.*cos(tt); height; ExParam*tt.*sin(tt)]';
    case 2  % Swiss Hole
        % Swiss Roll w/ hole example taken from Donoho & Grimes
        tt = (3*pi/2)*(1+2*rand(1,2*N));  
        height = 21*rand(1,2*N);
        kl = zeros(1,2*N);
        for ii = 1:2*N
            if ( (tt(ii) > 9) && (tt(ii) < 12))
                if ((height(ii) > 9) && (height(ii) <14))
                    kl(ii) = 1;
                end;
            end;
        end;
        kkz = find(kl==0);
        tt = tt(kkz(1:N));
        height = height(kkz(1:N));
        Samples = [tt.*cos(tt); height; ExParam*tt.*sin(tt)]';     
    case 3  % Corner Planes
        MyParam=1000;
        xMax = floor(sqrt(MyParam));
        yMax = ceil(MyParam/xMax);
        cornerPoint = floor(yMax/2);
        x=xMax*rand(1,N);
        y=yMax*rand(1,N);
        for k=1:N
            if y(k) <= cornerPoint
                X(k,:) = [x(k),y(k),0];
                ColorVector(k) = y(k);
            else
                X(k,:) = [x(k),cornerPoint+(y(k)-cornerPoint)*cos(pi*ExParam/180),(y(k)-cornerPoint)*sin(pi*ExParam/180)];
                ColorVector(k) = y(k);
            end;
            
        end
        Samples = X;
    case 4  % Punctured Sphere by Saul & Roweis
        xx=-5+10*rand(N,1);
        yy=-5+10*rand(N,1);
        rr2 = xx(:).^2 + yy(:).^2;
        [tmp ii] = sort(rr2);
        Y = [xx(ii(1:N))'; yy(ii(1:N))'];
        a = 4./(4+sum(Y.^2));
        Samples = [a.*Y(1,:); a.*Y(2,:); ExParam*2*(1-a)]';
    case 5  % Twin Peaks by Saul & Roweis
        xy = 1-2*rand(2,N);
        Samples = [xy; sin(pi*xy(1,:)).*tanh(3*xy(2,:))]';
        Samples(:,3) = ExParam * Samples(:,3);
    case 6  % Toroidal Helix by Coifman & Lafon
        t = (1:N)'/N;
        t = t.^(ExParam)*2*pi;
        Samples = [(2+cos(8*t)).*cos(t) (2+cos(8*t)).*sin(t) sin(8*t)];
    case 7  % Ueda's spiral
        t = (((0:N-1)/(N-1))*4*pi)';
        Samples = [(13-(0.5*t)).*cos(t) -(13-(0.5*t)).*sin(t) t];        
end

Samples=Samples';
Samples=Samples+StdNoise*randn(size(Samples));
