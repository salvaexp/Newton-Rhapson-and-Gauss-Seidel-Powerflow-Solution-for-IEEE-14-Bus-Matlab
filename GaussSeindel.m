%Algoritmo  para Gauss-Seidell
%Se utiliza para todo la ecuacion general del flujo de potencia.
%La tension de cada barra se calcula con la formula despejada.

%% Tension de cada barra en la iteracion actual, se obtiene el siguiente valor

 %Vi debe ser un vector de n entradas, n= # barras. los vectores son
 %horizontales
%clear all 
n=3; %# de barras
V= zeros(1,n); %Crea vector vacio de 1 fila, n columnas. Va a guardar todas las tensiones
V= V+1; %da por default el valor de cada tension=1

%Definimos los valores especiales que en este caso son V1 y V3
V(1)=1.05; 
V(3)=1.04;

Y= zeros(n,n); %se define Ybarra, matriz cuadrada. Aqui hay otros 100 pesos de programacion
%Datos iniciales:
Y= [20-50i -10+20i -10+30i ; -10+20i 26-52i -16+32i ; -10+30i -16+32i 26-62i];

%Datos que nos dan de potencia de las barras. Las desconocidas se dejarán en cero
P=[0 -4 2];
Q=[0 -2.5 0];

for h = 1:7 %%Conteo de las iteraciones, se va de 1 a 2
    
for e =2:n %se usará la letra "e" para representar la barra actual, como la letra i vista en el curso
    
    SUMi=0; %sumatoria en la tension
    for s = 1:n
        if (s==e)
            SUMi=SUMi; %no pasa nada
        else 
            SUMi= SUMi + Y(e,s)*V(s);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUMtot=0; %Suma total, se usa en calculo de Q
    for s = 1:n 
        SUMtot= SUMtot + Y(e,s)*V(s);
    end
    if e==3 %condicion en la que hay que calcular primero Q3
        Q(3)= -imag(conj(V(3))*SUMtot );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ya calculado Q3 se pasa a V3%%%%%%%
    
    V(e)=(1/Y(e,e))*((P(e)-Q(e)*1i)/(conj(V(e))) - (SUMi)); %Aqui se calculan las tensiones
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Reajuste del valor de V3=1,04
   if abs(V(3))== 1.04
       V(3)=V(3);
   else
       V(3)= sqrt(1.04^2 - imag(V(3))^2) + imag(V(3))*1i;
   end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end




