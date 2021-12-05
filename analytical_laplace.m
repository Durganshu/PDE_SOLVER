
%Evaluating the analytical expression for 2D Laplace Equation for temperature(x,y) =1 at
%right end and rest edges at 0.

x=linspace(0,100,100);  %Creates a vector between 0,1 with 100 divisions in between them
y=linspace(0,100,100);

H=100; L=100;

[X,Y] = meshgrid(x,y);

temperature=zeros(length(X),length(Y)); %Initialize the temperature

for n=1:200  
    coeff=(2*(1-cos(n*pi)))/((n*pi)*sinh((n*pi*L)/H));
    temperature_temp= (coeff)*(sinh(((n*pi)/H)*X)).*(sin(((n*pi)/H)*Y));
    temperature=temperature+temperature_temp;
end

 imagesc(x,y,temperature);
 colorbar;
 drawnow;

%contourf(X,Y,temperature);
%hold on;
