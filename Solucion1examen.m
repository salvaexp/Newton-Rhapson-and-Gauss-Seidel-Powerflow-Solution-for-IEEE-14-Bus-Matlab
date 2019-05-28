%///////////Primer Examen Sistemas de Potencia I \\\\\\\\\\\\\\\\
%Estudiante: Pablo Salvador Rodríguez Gómez, carné B45900
%Datos de la hoja en formato cdf:y
%Pot Base= 100.0 MVA
%Caso: IEEE 14 Bus Test Case


n=14; %# de barras
Y= zeros(n,n); %matriz vacía de admitancia de barras
Sbase=100;

%% Leer la tabla IMPORTANTE LEER ANTES DE EJECUTAR:
%Para leer los datos de la tabla de la ieee, esta se debe encontrar en la
%"carpeta actual" de Matlab, o debe ser indicada su ubicación completa 
%o relativa (su ruta o "path") en las variables "filename".
%El archivo de texto se llama 'ieee14cdf.txt' 

%Importar datos como tabla: tablaBus y tablaBranch:
%//////////////////primera tabla\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
filename = 'ieee14cdf.txt';
delimiter = ' ';
startRow = 3;
endRow = 16;
formatSpec = '%f%*s%*s%*s%*s%*s%f%f%f%f%f%f%f%*s%*s%f%f%f%f%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
tablaBus = table(dataArray{1:end-1}, 'VariableNames', {'Bus','Type','FinalVolt','FinalAngle','LoadMW','LoadMVAR','GenMW','GenMVAR','MaxMVAR','MinMVAR','ShuntCondG','ShuntSuscB'});
clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;

%/////////////////////////////segunda tabla\\\\\\\\\\\\\\\\\\\\\\\
filename = 'ieee14cdf.txt';
delimiter = ' ';
startRow = 19;
endRow = 38;
formatSpec = '%f%f%*s%*s%*s%*s%f%f%f%*s%*s%*s%*s%*s%f%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
tablaBranch = table(dataArray{1:end-1}, 'VariableNames', {'BarraOrigen','BarraDestino','Impedancia','Reactancia','Susceptancia','Vueltas'});
clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;

%\\\\\\\\\\\\\\\\\\\Fin creacion de tablas\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


%% Creación de la matriz Ybarra:
for h =1:20 %20 es el numero de filas en tablaBranch, es decir el numero de ramas.
    BarOr=tablaBranch{h,1};
    BarDes=tablaBranch{h,2};
    
    Res=tablaBranch{h,3};
    Reac=tablaBranch{h,4};
    Susc=tablaBranch{h,5};
    %se debe sumar la Susceptancia al inverso de la R+jX, G+jB
    Adm=(Res+1i*Reac)^(-1) + 1i*(Susc);
    
    Y(BarOr,BarDes)= -(Adm);
    Y(BarDes,BarOr)= -(Adm);
    

end
 %El cálculo de la admitacia propia de cada barra es: (Yparalelo + Yserie)
 %Yserie lo obtengo de la variable "Adm" o de las entradas xy
 %Yparalelo lo obtego de la tabla tablaBus

for g = 1:14 
        CondBus= tablaBus{g,11};
        SuscBus= tablaBus{g,12};
        
        SumaSerie=0;
        for f=1:14
            SumaSerie= SumaSerie + -1*(Y(g,f));
        end
        
        Y(g,g)= CondBus + 1i*SuscBus + SumaSerie;
end
%Matriz Y está lista



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cálculo de Gauss-Seidel
%Primero se va a guardar los datos iniciales de tensiones, angulos y potencias que nos dan.

Vini= zeros(1,n); 
for a=1:n %leo los datos iniciales de tensión dados, vector de 14.
    Vini(a)= tablaBus{a,3};
end
Angini= zeros(1,n); 
for a= 1:n % igual para ángulo inicial
    Angini(a)= tablaBus{a,4};
end

Angini = deg2rad(Angini); %se convierte de grados a radianes

V= zeros(1,n); %Esta variable va a guardar las tensiones complejas de cada bus. Parte con el valor inicial.
for a= 1:n
    V(a)= Vini(a).*exp(1i*Angini(a));
end


%Se definen potencias iniciales: suma de lo generado más lo consumido
P= zeros(1,n); 
for a=1:n
    P(a) = tablaBus{a,7} - tablaBus{a,5};
end
Q = zeros(1,n);
for a=1:n
    Q(a) = tablaBus{a,8} - tablaBus{a,6};
end

Q=Q/Sbase;

P=P/Sbase;

%Segun la teoria, hay 3 tipos de barra:
%Oscilante[no se hacen calculos de tension] tipo 3
%carga(P-Q)[se calcula V y theta] tipo 0
%generación(P-V) [se calcula theta y Q] (el límite causa que se transforme
%en barra P-Q) tipo 2
BusType= zeros(1,n);
for a= 1:n
    BusType(a)= tablaBus{a,2}; %En este vector se va guardando el tipo de barra
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Primera parte: cálculo de las 14 tensiones
iteracion= 1;
fin=2; %se usan estas dos variables para controlar el error y las iteraciones.
while(iteracion<fin) %comienzo de iteraciones. Si se llega al error deseado, se detiene el loop.
    
    %funcion de error: si el error no es menor a 10^-4, las iteraciones siguen
    CuentaError= V; %se verá el error en todas las barras. En esta variable se guarda el valor actual.
    
    %Se procede al cálculo de tensiones: según la fórmula dada en el curso
    for k= 2:n %k=barra actual. se comienza desde 2 porque 1 es barra osc.
        
        SUMk=0; %sumatoria para la formula de cálculo de tension de cada barra
        for s=1:n
            if  (s==k)
                %no hacer nada
            else 
                SUMk=SUMk + Y(k,s)*V(s);
            end
        end %fin sumatoria

        %% Si la barra es P-V: no tenemos Q y se debe cumplir valor de V dado
        
        if (BusType(k)==2)
            SUMtot=0; %variable para guardar la suma de corrientes "Y*V"
            for s= 1:n
              SUMtot = SUMtot + Y(k,s)*V(s);
            end 
        Q(k)= -imag(conj(V(k))*SUMtot); %obtiene Q de la barra "k"
        end
        %% Obtención de V
        %fórmula de flujo de potencia despejada para tensión de barra:
        V(k)= (1/Y(k,k))*( (P(k)-1i*Q(k))/(conj(V(k))) - (SUMk) ); 
        
        %Si la barra es P-V, se debe reajustar al V dado:
        if (BusType(k)==2)
            if (abs(V(k))~= abs(Vini(k)))
               V(k) = sqrt(  Vini(k)^2 - imag(V(k))^2 ) + imag(V(k))*1i; 
            end   
        end %if barra PV
    end  %for de barras 1 a 14
    
    
   error= abs(V-CuentaError);
   for i= 1:n
       if ((error(i))>(10^(-4)))
           fin=fin+1;
           break
       end
   end
   %el while continua si iteracion<fin. el valor de "iteracion" siempre aumenta, 
   %el valor de "fin" sólo cuando el error no es suficientemente bajo.
    iteracion=iteracion+1;
end %iteraciones
disp('El numero de iteraciones para el metodo G-S fue:');
display(iteracion-1);
%disp('Las tensiones resultantes de cada barra son');
%disp(V);
 
%%%%%%%%%%%%%%%%%%%%%%%FIN DEL CÁLCULO DE TENSIONES%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Segunda parte: cálculo de las Potencias P1 y Q1 generadas por barra 1
%Se utiliza la fórmula general del flujo de potencia, ya con las tensiones
%conocidas: Si*= (Vi*).*(SumYikVk)

%De las barras de tipo 2 y 3 ya se tiene el valor de S.
SUMtot=0; %variable para guardar la suma de corrientes "Y*V"
for k= 1:n
   SUMtot = SUMtot + Y(1,k)*V(k);
end 
    
Sconj= conj(V(1))*SUMtot;
P(1)= real(Sconj);
Q(1)= -imag(Sconj);
%disp('La barra generadora 1 produce MVA:');
%display((P(1)+1i*Q(1))*100); %el display de todas las barras se hace despues





%Código para avisar cuando un generador alcance su límite de potencia: aplica para
%barras de tipo P-V

MaxMVAR= zeros(1,n);
for i= 1:n
    MaxMVAR(i)= tablaBus{i,9};
end

MinMVAR= zeros(1,n);
for i= 1:n
    MinMVAR(i)= tablaBus{i,10};
end


for i= 1:n
    if (BusType(i)==2)
        if (  (Q(i)-(tablaBus{i,6}/100))  >(MaxMVAR/100)) %limite max de potencia es menor a potencia generada
            disp('El generador ha alcanzado su límite máximo de potencia reactiva');
        end   
            
        if(  (Q(i)-(tablaBus{i,6}/100))  <(MinMVAR/100))  %no se produce suficiente
            disp('El generador ha alcanzado su límite mínimo de potencia reactiva');
        end
        
%         if(Q(i)>MinMVAR) 
%             disp('El generador NO ha alcanzado su límite mínimo de potencia reactiva');
%             disp('actualmente se genera:');
%             disp(Q(i));
%         end
        
        %disp('The result is:')
        %disp(P(i))
    end
end

%%%%%%%%%%%%%%%%%%FIN DE CALCULO DE POTENCIAS DE BUS 1%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Tercera parte: usando tensiones de cada barra, cálculo de las corrientes 
%de línea y potencias de línea.

%Corrientes a través de las líneas:
%De barra i a barra j: Iij= Iserie+Ipar= Ys(Vi-Vj)+Ypi*Vi
%Ys: datos ij. (informacion de ramas)
%Yp: datos ii. (informacion de buses) 


Yp=zeros(1,n);
for g = 1:14 
        CondBus= tablaBus{g,11};
        SuscBus= tablaBus{g,12};
        Yp(1,g)= CondBus + 1i*SuscBus; %se guarda información de impedancia de derivación de cada bus
end

%%%%Se comenzará el cálculo de corrientes:
I=zeros(n,n);
Ibruto=zeros(n,n);
%el resultado de las corrientes será  una matriz no simétrica con diagonal 0

for i=1:n %Para la corriente de la barra i...
    for j=1:n %...hacia la barra j
        Ibruto(i,j)= -Y(i,j)*(V(i)-V(j))+Yp(i)*V(i);
    end
    Ibruto(i,i)=0; %Se hace esto porque no tenemos definida la corriented e una barra a sí misma.
end

%Se tiene una matriz "Ibruto" que guarda corrientes entre todas las barras.
%Esto es incorrecto pues no todas están conectadas físicamente entre sí.
%%%%%La Correccion a esto: matriz "I"
%Se dejarán intectas las corrientes entre barras definidas en la tabla tablaBranch y las demás se borran. 

for h =1:20 %20 es el numero de ramas
    BarOr=tablaBranch{h,1};
    BarDes=tablaBranch{h,2};
    
    I(BarOr,BarDes)=Ibruto(BarOr,BarDes);
    I(BarDes,BarOr)=Ibruto(BarDes,BarOr);
end

%El resultado de las corrientes se guada en matriz "I" y se lee: Corriente
%desde barra "fila" hacia barra "columna".


%///////////////////////////////////Flujo de potencia\\\\\\\\\\\\\\\\\\\\\
%Flujo de potencia: Potencia que envía cada barra a las demas:
%Se calcula como: Sij=Vi*conj[Iij]

FlujoS=zeros(n,n);  %Es la tensión de la barra i, por la corriente que sale de ella i->j
for i=1:n %Flujo desde barra i...
    for j=1:n %...hasta barra j
        FlujoS(i,j)= V(i)*conj(I(i,j));
    end
end
 % El resultado de los flujos de potencia se guarda en la matrizz FlujoS, y
 % se lee: El flujo desde la barra "fila" hacia la barra "columna"
 
 
%////////////////////////Pérdidas en lineas\\\\\\\\\\\\\\\\\\\\\\
% Pérdidas en las líneas: la diferencia del flujo de potencia en los
% extremos de cada rama.

PerdidasS=zeros(n,n); %Se suma el flujo de potencia de las dos barras en una rama.
for i=1:n
    for j=1:n
        PerdidasS(i,j)=FlujoS(i,j)+FlujoS(j,i);
    end
end

%Se obtienen las pérdidas en matriz "PerdidasS". Es una matriz simétrica, se lee: pérdida entre
%barra "fila" y barra "columna", y viceversa.

%/////////////////////////FIN CÁLCULO POTENCIAS\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\







%% Lectura del estado de las barras:
%Esta seccion es para imprimir en pantalla los resultados anteriores.

for i=1:n
    fprintf('La potencia total en la barra %d es: \n', i);
    disp((P(i)+1i*Q(i))*100);
    %disp('La potencia total en la barra ', i, ' es ',  )
end

for i=1:n
    if(BusType(1)==2 || BusType(1)==3)
        fprintf('La potencia generada en la barra %d es: \n',i);
        disp( (P(i)+1i*Q(i))*100 + tablaBus{i,5}+ 1i*tablaBus{i,6} );
        %disp('La potencia generada en la barra ', i, ' es ', (P(i)+1i*Q(i))*100 + tablaBus{1,5}+ 1i*tablaBus{1,6})
    end
end

for i=1:n
    fprintf('Según el método de G-S, la tensión en la barra %d es %f V con angulo de %f °\n',i , abs(V(i)) , rad2deg(angle(V(i)))  );
    %disp('La tensión en la barra ', i, ' es ', abs(V(i)) ,' con angulo de ', rad2deg(angle(V(i))) )
end


for h =1:20 %20 es el numero de ramas
    BarOr=tablaBranch{h,1};
    BarDes=tablaBranch{h,2};
    
    fprintf( 'El flujo de potencia de la barra %d a la barra %d , es de: \n', BarOr , BarDes );
    disp((FlujoS(BarOr,BarDes)*100));
    %disp('El flujo de potencia de potencia de la barra ', BarOr,' a la barra ',BarDes,' es de ', (FlujoS(BarOr,BarDes)*100));
    
    fprintf( 'El flujo de potencia de la barra %d a la barra %d , es de: \n', BarDes , BarOr );
    disp((FlujoS(BarDes,BarOr)*100));
    
end

for h =1:20 %20 es el numero de ramas
    BarOr=tablaBranch{h,1};
    BarDes=tablaBranch{h,2};
    
    fprintf('Las pérdidas en la rama %d a %d son de: \n',BarOr,BarDes );
    disp( (PerdidasS(BarOr,BarDes)*100) );
    %disp('Las pérdidas en la rama desde ', BarOr,' a la barra ',BarDes,' son de ', (PerdidasS(BarOr,BarDes)*100));
    
    fprintf('Las pérdidas en la rama %d a %d son de: \n',BarDes,BarOr );
    disp( (PerdidasS(BarDes,BarOr)*100) );
end






%% Cálculo de Newton-Rhapson
%Este método utiliza variables de estado (anglulo, tension), variables de
%control (P y Q) y un jacobiano, que se obtiene de derivar formulas de potencia
%respecto a cada variable de estado.

%El procedimiento de este método consiste en tres pasos por iteración:
%Primero, obtener vector de control y obtener el error de control con datos reales, luego obtener el Jacobiano 
%Segundo se obtiene el error de las variables X a partir del Jacob invertido junto al error de control
%Tercero se obtiene el valor final de X a partir del valor inicial y el error

%Dependiendo del tipo de barra, las variables de estado son diferentes.
%Para barras P-V(tipo 2) la variable es angulo (hay que obtener al final la
%potencia Q) y control U es con P
%Para barras P-Q(tipo 0) la variable es angulo y tension, y control U es con
%P y Q.

%Algoritmo: para la iteracion actual 
%errorU= Uini - U(generado con f y Xanterior.)
%errorX= inv[J]*errorU
%X = Xanterior + errorX

%Hay que tener fórmulas explícitas para cada variable de control. "f"
%Hay que obtener el Jacobiano con estas fórmulas. "J"

%Datos a reutilizar: Y, Vini, Angini, BusType
%Datos a generar de nuevo: Potencias, Tension a calcular.
clearvars P Q V;

Pini= zeros(1,n); 
for a=1:n
    Pini(a) = tablaBus{a,7} - tablaBus{a,5};
end
Qini = zeros(1,n);
for a=1:n
    Qini(a) = tablaBus{a,8} - tablaBus{a,6};
end
Qini=Qini/Sbase;
Pini=Pini/Sbase;
%%%%%%%%%%%%%%%%%%
%Obtener propiedades de la matriz Y
angY=angle(Y);
magY=abs(Y);
%%%%%%%%%%%%%%%%%%%

%Ordenamiento de variables: X=primero angulos y luego tensiones/ U=primero P y luego Q /


%Para generar un U cada iteración se necesitan las funciones:
%Pi= SUMk Yik*Vk*Vi*cos(angk-angi+angYik)
%Qi= -SUMk Yik*Vi*Vk*sen(angk-angi+angYik)
%Estas se crearan a continuacion:


% Se crean vectores de variables de tension y fase, para tener formulas de potencia:
V = sym('v', [1 n]); %crea un vector simbólico al cual se puede acceder a sus entradas
ang = sym('ang', [1 n]);

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%En esta seccion se van a guardar angulos, potencia real, y fórmulas de potencia real
%Este codigo se corre una sola vez.
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
tamano=0;
for i=1:n %este "for" es para obtener el tamaño a usar.
    if (BusType(i)==0) || (BusType(i)==2)
        tamano=tamano+1;
    end
end 


angEva= zeros(tamano,1); %para guardar el valor inicial de X, para la primera iteracion
angX= sym('angx',[1 tamano]) ;%para guardar las variables de angulo que se usaran para derivar.

P= zeros(1,tamano);
fp=sym('fp', [1 tamano]); %crea vector para guardar variables simbolicas, seran las funciones de fp


identiBarraPANG=zeros(1,tamano); % Para tener guardado a qué barra pertenecen angulos y potencias reales


posicion=0;
for i=1:n
    if (BusType(i)==0) || (BusType(i)==2)
        posicion=posicion+1;
        angEva(posicion)=Angini(i); %angEva guarda el valor inicial de las variables de estado de angulo
        P(posicion)=Pini(i);
        
        identiBarraPANG(posicion)=i; 
        
        %variables para jacobiano
        angX(posicion)=ang(i);
        
        %generar funciones fp:
        SUMA=0;
        for k=1:n
            SUMA= SUMA+ magY(i,k)*V(i)*V(k)*cos(ang(k)-ang(i)+angY(i,k));
        end
        fp(posicion)=SUMA;
        
    end
end

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%guardar tensiones, Q y formulas de potencia Q:
%Este codigo tambien solo se corre una vez.
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
tamano=0;
for i=1:n %Otro for para obtener el tamaño a usar
    if (BusType(i)==0)
        tamano=tamano+1;
    end
end

Veva= zeros(tamano,1); %Veva es el valor inicial que tienen las tensiones de las barras P-V (ultimas variavles de X)
VX= sym('Vx',[1 tamano]);%para guardar las variables de tension que se usaran para derivar.

Q=zeros(1,tamano); 
fq=sym('fq', [1 tamano]); %crea vector para guardar funciones


identiBarraQV=zeros(1,tamano); % para tener guardado a qué barra pertenece este valor de Veva y de Q


posicion=0;
for i=1:n
    if (BusType(i)==0)
        posicion=posicion+1;
        Veva(posicion)=Vini(i); %se generan los valores iniciales de X
        Q(posicion)=Qini(i);
        
        identiBarraQV(posicion)=i;
        
        
        %variables para jacobiano
        VX(posicion)=V(i);
        
        %generar funciones fq:
        SUMA=0;
        for k=1:n
            SUMA= SUMA+ magY(i,k)*V(i)*V(k)*sin(ang(k)-ang(i)+angY(i,k));
        end
        fq(posicion)=-SUMA;
    end
end





%///////\\\\\\\\\////  Aqui se concatena lo generado anteriormente  \\/////////////\\\\\\\\\///////
%En el caso de las 14 barras: Vectores de U y X son tamaño 13+9=22
Xini= [angEva;Veva];        %Aqui guadramos valores de X que seusan al final de cada iteracion    
Uini= [P,Q];                %Valores de U que se usaran de referencia en todas las iteraciones.
f=[fp,fq];                  %funciones simbólicas
X=[angX,VX]; % Importante: aqui guarda cuales variables usar para derivar. Se usa en Jacobiano.

%///////\\\\\\\\\\//////////////\\\\\\\\\\//////////////\\\\\\\\\\//////////////\\\\\\\\\\///////
%Ya tengo X y U iniciales

%% Creacion del jacobiano: 

%Jacobiano: resultado de derivar parcialmente una función por fila, una 
%variable por columna, y luego evaluar [V,ang] en [Vini,Angini]:

[~,columnas]= size(Uini); %obtengo el tamaño que llevará el jacobiano
J=sym('J',[columnas columnas]); %será matriz simbólica, como f.

filas=columnas;

for j=1:filas %para recorrer las filas
    for k= 1:columnas %para recorrer las columnas
        J(j,k)= diff(f(j),X(k));%cada entrada será la derivada de f(j) respecto a X(k)
    end
end
%\\\\\\\\\J se demora un rato en calcular, pero solo hay que hacerlo una vez en todo la simulacion/////////////

%% Comienzo de iteraciones:




Viter= Vini; %Estos dos vectores llevaran los datos actualizados de tensiones y angulos
Angiter= Angini;

iteracion= 1;
fin=2; %se usan de nuevo estas dos variables para controlar el error y las iteraciones.

while (iteracion < fin)




%% primer paso de iteracion: obtencion de errorU= Uini - U(generado con f y Xanterior)
%Tengo Xini, que guardael valor de iteracion anterior. Uini guarda Potencias
%f guarda las formulas. 

%sintaxis de matlab: si hago subs(f), obtengo el resultado de evaluar todas las variables que hay en f.
%f sigue guardando las expresiones de potencia. ejemplo: Ugenerado=subs(f)
%Para esto debo darle valores a todas las variables en los vectores simbolicos: V() y ang().

errorU= Uini - subs(f,[V,ang],[Viter,Angiter]); %la funcion subs evalua sobre f, reemplazar V por "V iter" y ang por "Ang iter".
errorUvpa= vpa(errorU); %para pasar a formato double.

%% segundo paso: obtencion de errorX= inv[J]*errorU
%Se usa el jacobiano ya calculado. Aqui se evalua con V iter, Ang iter.

Jo=subs(J,[V,ang],[Viter,Angiter]); %CALCULA JACOBIANO en cada iteracion
Jvpa=vpa(Jo); %se pasa Jacobiano a formato double.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Al hacer cálculos con matrices symbolics, Matlab se vuelve loco. Se
%pasará a matrices normales. Se llamaran "XXXmat"
[filas,columnas]=size(errorUvpa);
errorUmat=zeros(columnas, filas);
for j=1:filas %crea las columnas j
    for i=1:columnas %crea las filas i
        errorUmat(i,j)=errorUvpa(j,i);
    end
end

[filas,columnas]=size(Jvpa); %se pasa de un arreglo sym a matriz.
Jmat=zeros(filas,columnas);
for i=1:filas
    for j=1:columnas
        Jmat(i,j)=Jvpa(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


errorX=Jmat\errorUmat; %Jmat\errorUmat es igual que inv(Jmat)*errorUmat


%% tercer paso: obtención de X = Xanterior + errorX
%En este paso se actualiza el valor de X que será tomado en la siguiente
%iteración.

%Los vectores Viter y Angiter guardan los valores de Tension y ang de cada iteracion.  
%Xini va a recibir el resultado de la simulación, y lo va a llevar al nuevo valor de Viter y Angiter. 

Xini=Xini + errorX;

%El dato de esta iteración que se usará en errorU(f) y en J para la siguiente,
%se guarda en el vector Viter, Angiter
    
    
    
    [tamanoANG,~]= size(angEva);
    [tamanoV,~]= size(Veva);
   	
    for i=1:tamanoANG %del 1 al 13
        indexANG= identiBarraPANG(i); %guarda la barra que corresponde a la entrada i de AngEva
        Angiter(indexANG)= Xini(i);
    end

    for i= 1:tamanoV %del 1 al 9
        indexV= identiBarraQV(i); %guarda la barra que corresponde a la entrada i de Veva
        Viter(indexV)= Xini(i+ tamanoANG); %continua leyendo Xini donde lo había terminado.
    end

%%%%%%%%%%%%Cuenta del error para hacer iteraciones%%%%%%%%%%%%%%%%%%%
    for i= 1:columnas
       if (abs(errorX(i))>(10^(-4)))
           fin=fin+1;
           break
       end
    end
   %el while continua si iteracion<fin. el valor de "iteracion" siempre aumenta, 
   %el valor de "fin" sólo cuando el error no es suficientemente bajo.
    iteracion=iteracion+1;
%%%%%%%%%%FIN%%%%%%%%%%%%

end

Vnr= zeros(1,n); %Esta variable va a guardar las tensiones complejas de cada bus.
for a= 1:n
    Vnr(a)= Viter(a).*exp(1i*Angiter(a));
end



disp('El numero de iteraciones para el método N-R fue:');
disp(iteracion-1);
%disp('Las tensiones resultantes son:');
%disp(Vnr);

for i=1:n
    fprintf('Según el método de N-R, la tensión en la barra %d es %f V con angulo  de %f ° \n',i , Viter(i) , rad2deg(Angiter(i))  );
    %disp('La tensión en la barra ', i, ' es ', abs(V(i)) ,' con angulo de ', rad2deg(angle(V(i))) )
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Seccion adicionada: potencias calculadas a partir de las tensiones obtenidas con N-R: copy paste de G-S
%SUPONGO (porque en clase nunca lo explicaron) que lo correcto para esta parte es
%obtener las potencias a partir de las ecuaciones de sumatoria
%Yik*Vk*Vi*cos(dk-di+yik). Pero para lograrlo habría que modificar mucho
%este código, y tengo que entregar esto en 9 horas al profe.

clearvars P Q V Sconj I Ibruto FlujoS PerdidasS;



%Se definen potencias iniciales DE NUEVO (G-S)
P= zeros(1,n); 
for a=1:n
    P(a) = tablaBus{a,7} - tablaBus{a,5};
end
Q = zeros(1,n);
for a=1:n
    Q(a) = tablaBus{a,8} - tablaBus{a,6};
end

Q=Q/Sbase;

P=P/Sbase;


%Segunda parte: cálculo de las Potencias P1 y Q1 generadas por barra 1
%Se utiliza la fórmula general del flujo de potencia, ya con las tensiones
%conocidas: Si*= (Vi*).*(SumYikVk)

%De las barras de tipo 2 y 3 ya se tiene el valor de S.
SUMtot=0; %variable para guardar la suma de corrientes "Y*V"
for k= 1:n
   SUMtot = SUMtot + Y(1,k)*Vnr(k);  %%A partir de aqui la tension Vnr es la obtenida del método N-R
end 
    
Sconj= conj(Vnr(1))*SUMtot;
P(1)= real(Sconj);
Q(1)= -imag(Sconj);
%disp('La barra generadora 1 produce MVA:');
%display((P(1)+1i*Q(1))*100); %el display de todas las barras se hace despues










% Tercera parte: usando tensiones de cada barra, cálculo de las corrientes 
%de línea y potencias de línea.

%Corrientes a través de las líneas:
%De barra i a barra j: Iij= Iserie+Ipar= Ys(Vi-Vj)+Ypi*Vi
%Ys: datos ij. (informacion de ramas)
%Yp: datos ii. (informacion de buses) 


Yp=zeros(1,n);
for g = 1:14 
        CondBus= tablaBus{g,11};
        SuscBus= tablaBus{g,12};
        Yp(1,g)= CondBus + 1i*SuscBus; %se guarda información de impedancia de derivación de cada bus
end

%%%%Se comenzará el cálculo de corrientes:
I=zeros(n,n);
Ibruto=zeros(n,n);
%el resultado de las corrientes será  una matriz no simétrica con diagonal 0

for i=1:n %Para la corriente de la barra i...
    for j=1:n %...hacia la barra j
        Ibruto(i,j)= -Y(i,j)*(Vnr(i)-Vnr(j))+Yp(i)*Vnr(i);
    end
    Ibruto(i,i)=0; %Se hace esto porque no tenemos definida la corriented e una barra a sí misma.
end

%Se tiene una matriz "Ibruto" que guarda corrientes entre todas las barras.
%Esto es incorrecto pues no todas están conectadas físicamente entre sí.
%%%%%La Correccion a esto: matriz "I"
%Se dejarán intectas las corrientes entre barras definidas en la tabla tablaBranch y las demás se borran. 

for h =1:20 %20 es el numero de ramas
    BarOr=tablaBranch{h,1};
    BarDes=tablaBranch{h,2};
    
    I(BarOr,BarDes)=Ibruto(BarOr,BarDes);
    I(BarDes,BarOr)=Ibruto(BarDes,BarOr);
end

%El resultado de las corrientes se guada en matriz "I" y se lee: Corriente
%desde barra "fila" hacia barra "columna".


%///////////////////////////////////Flujo de potencia\\\\\\\\\\\\\\\\\\\\\
%Flujo de potencia: Potencia que envía cada barra a las demas:
%Se calcula como: Sij=Vi*conj[Iij]

FlujoS=zeros(n,n);  %Es la tensión de la barra i, por la corriente que sale de ella i->j
for i=1:n %Flujo desde barra i...
    for j=1:n %...hasta barra j
        FlujoS(i,j)= Vnr(i)*conj(I(i,j));
    end
end
 % El resultado de los flujos de potencia se guarda en la matrizz FlujoS, y
 % se lee: El flujo desde la barra "fila" hacia la barra "columna"
 
 
%////////////////////////Pérdidas en lineas\\\\\\\\\\\\\\\\\\\\\\
% Pérdidas en las líneas: la diferencia del flujo de potencia en los
% extremos de cada rama.

PerdidasS=zeros(n,n); %Se suma el flujo de potencia de las dos barras en una rama.
for i=1:n
    for j=1:n
        PerdidasS(i,j)=FlujoS(i,j)+FlujoS(j,i);
    end
end

%Se obtienen las pérdidas en matriz "PerdidasS". Es una matriz simétrica, se lee: pérdida entre
%barra "fila" y barra "columna", y viceversa.

%/////////////////////////FIN CÁLCULO POTENCIAS\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



%% Lectura del estado de las barras:
%Esta seccion es para imprimir en pantalla los resultados anteriores.

fprintf('RESULTADOS DE NEWTON-RHAPSON');

for i=1:n
    fprintf('La potencia total en la barra %d es: \n', i);
    disp((P(i)+1i*Q(i))*100);
    %disp('La potencia total en la barra ', i, ' es ',  )
end

for i=1:n
    if(BusType(1)==2 || BusType(1)==3)
        fprintf('La potencia generada en la barra %d es: \n',i);
        disp( (P(i)+1i*Q(i))*100 + tablaBus{i,5}+ 1i*tablaBus{i,6} ); %% CORREGIR A tablaBus{i,7}+ 1i*tablaBus{i,8}
        %disp('La potencia generada en la barra ', i, ' es ', (P(i)+1i*Q(i))*100 + tablaBus{1,5}+ 1i*tablaBus{1,6})
    end
end

% Esto ya fue mostrado antes, da exactamente igual.
% for i=1:n
%     fprintf('Según el método de N-R, la tensión en la barra %d es %f V con angulo de %f °\n',i , Viter(i) , rad2deg(Angiter(i))  );
%     %disp('La tensión en la barra ', i, ' es ', abs(Vnr(i)) ,' con angulo de ', rad2deg(angle(Vnr(i))) )
% end


for h =1:20 %20 es el numero de ramas
    BarOr=tablaBranch{h,1};
    BarDes=tablaBranch{h,2};
    
    fprintf( 'El flujo de potencia de la barra %d a la barra %d , es de: \n', BarOr , BarDes );
    disp((FlujoS(BarOr,BarDes)*100));
    %disp('El flujo de potencia de potencia de la barra ', BarOr,' a la barra ',BarDes,' es de ', (FlujoS(BarOr,BarDes)*100));
    
    fprintf( 'El flujo de potencia de la barra %d a la barra %d , es de: \n', BarDes , BarOr );
    disp((FlujoS(BarDes,BarOr)*100));
    
end

for h =1:20 %20 es el numero de ramas
    BarOr=tablaBranch{h,1};
    BarDes=tablaBranch{h,2};
    
    fprintf('Las pérdidas en la rama %d a %d son de: \n',BarOr,BarDes );
    disp( (PerdidasS(BarOr,BarDes)*100) );
    %disp('Las pérdidas en la rama desde ', BarOr,' a la barra ',BarDes,' son de ', (PerdidasS(BarOr,BarDes)*100));
    
    fprintf('Las pérdidas en la rama %d a %d son de: \n',BarDes,BarOr );
    disp( (PerdidasS(BarDes,BarOr)*100) );
end




