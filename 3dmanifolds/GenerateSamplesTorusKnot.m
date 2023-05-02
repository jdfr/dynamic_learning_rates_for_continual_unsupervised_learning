function [Muestras]=GenerateSamplesTorusKnot(p,q,NumMuestras)

Phi=2*p*pi*rand(1,NumMuestras);

Muestras=zeros(3,NumMuestras);

Muestras(1,:)=(2+cos(q*Phi/p)).*cos(Phi);
Muestras(2,:)=(2+cos(q*Phi/p)).*sin(Phi);
Muestras(3,:)=sin(q*Phi/p);

Muestras=Muestras+0.01*randn(size(Muestras));
