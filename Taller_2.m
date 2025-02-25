% clc; clear; close all;


sigma=1e-5; % Espesor de la pelicula estancanda de B
C_A_0=1; %Concentración inicial
D_AB=5.5819e-10; % Difusividad [mol/m^2 s]
V=1; % Volumen total del líquido
n=300;

%% Valores usados en corridas de prueba
%k=linspace(8e-9, 2e1, 20); %Constante de velocidad de reacción
%C_A_front=[1, 0.5]; % Condiciones frontera [inical, final de la capa de reacción]


%  S_b=[1e0, 1e-2, 1e-3, 1e-4, 1e-5];
%     for j=1:5
%        S=S_b(j)

%% La superficie de la burbjas se varía manualmente para evitar los errores
%   del metodo númerico en la grafica. 
  
     S=1e-8; %Superficie de las burbujas
%% El valor de phi se varia manualmente para cada vector de superficie    
%  phi=linspace(1e-3, 1e-10, n);
phi=1e-6:1e-7:1e-5

% Vector de phi que se usó para encontrar los rangos de validos para el metodo numerico     
%  phi=[1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001 0.000005, 0.000001]
 
 C_A=zeros(length (phi),100); % Iniciación de matriz composición
 u=zeros(length (phi),100); % Iniciación de matriz dC_A/dr
 zmesh=linspace (0,sigma,5); % Mallado
        
% Se calcula la constante de velocidad de reacción homogenea según
% la definición de phi=sqrt(k.*sigma^2./D_AB).
  k=phi.^2.*D_AB./sigma^2; % Constante de velocidad de reacción [1/s]
        
for i=1:length (phi)
    
 % Según ecuación 18.4-15 Bird et al. 2ed.
 % B, la relación C_A_sigma/C_A_0
   B=1/(cosh(phi(i))+((V/S/sigma)*phi(i)*sinh(phi(i)))); 
    
  % Apartir de la relación B se calcula la concentración en el final de la
  % pelicula reactiva
    C_Asig=B*C_A_0;
       
  % Por lo que las condiciones frontera quedan
    C_A_front=[1, C_Asig];
    
  %% Solución del sistema de ecuaciones
  solinit = bvpinit (zmesh, @(z)guess(z,sigma,C_A_front,D_AB,k(i)));
    
   sol=bvp5c(@(z, Y)odefun_F(z,Y,k(i),D_AB, V, S),@(ya,yb)bcfun(ya,yb,C_A_front), solinit);

   z=linspace(0,sigma);
   Y=deval(z,sol);
   C_A(i,:)=Y(1,:); % Concentraciones
   u(i,:)= Y(2,:); % dC_A/dr
   N_mod(i)=(phi(i)/sinh(phi(i)))*(cosh(phi(i))-(1/(cosh(phi(i))+(V/S/sigma)*phi(i)*sinh(phi(i)))));
        
end

   %phi=sqrt(k.*sigma^2./D_AB);
    
%Calculo de la velocidad de absorción en r=0
    flux=u(:,2).*(-sigma/C_A_front(1));

%% Grafico de velocidad de absorción vs phi
loglog(phi, flux', 'b', phi, N_mod,'--b');
hold on

legend ('V/S\delta=','V/S\delta=')
title('perfil flux')
hold on

%    end
   



% figure(2)
% plot(z*100,C_A(1,:));
% hold on
% plot(z*100,C_A(7,:));
% plot(z*100,C_A(19,:));
% plot(z*100,C_A(20,:));
% title('perfil concentración vs velocidad de absorción')




