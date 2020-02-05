function [PotLibre,Pmwm] = Cost321MWM()
%Modelado de propagacion de onda de radiofreceuncia en el rango 2-5 GHz
%PotLibre es la potencia recibida en funcion del tranmisor dependiente de
%la distancia. Pmwm concidera paredes 
%Perdidas en espacio libre Lfs segun  UIT-R P.1238
%% 
%Caracteristicas del dispositivo transmisor y receptor  
Pt=30; %Potencia de trasmision dBm
Gt=3; %Ganancia de antena transmisora dBi
Pr=-75; %Ganancia de la antena receptora dBi
Gr=3; %Ganacia receptor 3dBi
f=2.4e9; %Frecuencia de transmision 
lambda=3e8/f;
resolucion=0.05;%RResolucion de calculo 5cm
dmax=10;%auxilar radio maxima de caculo 
%Contantes empiricas para modelado de interiores 
%Aragón-Zavala, A. (2017). Indoor Wireless Communications: From Theory to Implementation. John Wiley & Sons.
Lc=0.5; %Constantes de perdida guasiana  
Kw=1; %Numero de paredes 
n=0;%Numero de suelos 
Lw=3.4; %Perdida debido a muros delgado, sTabla
Lf=18.3;% Baldosa de hormigo 30cmm espesor 
Kf=0;% Numero de suelos 
b=0.46;% Parametro empirico 
%%
%Mapa y matriz de paredes 
%Lectura del archivo cad
%author de la funcion f_LectDxf = SebaTM 
 %Jul 27 2009
[c_Line,c_Poly,c_Cir,c_Arc,c_Poi] = f_LectDxf('Plano1.dxf');
 %Codigo propio, solo se trabajo con lineas [px,py] son los puntos de las coordenadas
 %tridimensionales de Autocad a matlab 
 [px py]=limites(c_Line);
 px=round(px,2);
 py=round(py,2);
%%
%Generar una matriz con los puntos de las paredes 
dmax=max(max(px),max(py));
dmin=min(min(px),min(py));
x=dmin:resolucion:dmax;
y=dmin:resolucion:dmax; 
matrizPlano=zeros(length(y),length(x));
xn=x;yn=y;
%Ubicar las paredes en la matriz matrizPlano. 
matrizPlano=ubicarPuntos(matrizPlano,x,flip(y),px,py);
%Ubicacion de transmisor y receptor 
[x,y] = getpts;
Tx=[round(x),round(y)];
figure(2)
hold on
plot(Tx(1),Tx(2),'c*');
figure(4) 
imshow(matrizPlano);
title('matriz puntos');
figure(5)
%%
%{
Calculos basados en: 
Indoor Wireless Communications
From Theory to Implementation
Alejandro Aragón-Zavala
5.5.9 COST-231 Multiwall Model pag 126 en adelante 
%}
% 1 Perdidas en espacio libre 
%  Malla de distancias , ubicacion del router Tx
lizq=resolucion:resolucion:(Tx(1)-dmin);
lup=resolucion:resolucion:(dmax-Tx(2));
lder=resolucion:resolucion:(dmax-Tx(1));
ldown=resolucion:resolucion:(Tx(2)-dmin);
Rango=[flip(lizq) 0 lder];
Domino=[flip(ldown) 0 lup];
for i=1:length(Rango)
    for j=1:length(Domino)
        m(j,i)=sqrt(Rango(i)^2+Domino(length(Domino)-j+1)^2);
    end 
end 
%1 Perdidas en espacio libre Lfs segun  UIT-R P.1238
%Para interiores 
N=30;%Coeficiente de perdidade potencia por distancia ITU 
Lfn=0;%Factor de perdidas por penetracion en funcion de numero de pisos
Lfs=20*log10(f/1e6)+N*log10(m)+Lfn-28;%28 en dB;
Lfs(~isfinite(Lfs))=0;
Prxlibre=Pt+Gt+Gr-Lfs;
%Grafica perdidas espacio libre lineal
figure(1);
x=resolucion:resolucion:dmax;
x=[-flip(lizq) 0 lder];
[n,r]=size(Prxlibre);
p=find(Prxlibre((n-1)/2,:)==max(Prxlibre((n-1)/2,:)));
p=find(Prxlibre(:,p)==max(Prxlibre(:,p)));
plot(x,Prxlibre(p,:));
title('Modelos de propagacion');
PotLibre=Prxlibre;
%%
% 2 Perdidas debido a paredes. matrizPlano y m
Lmuros=DeteccionMuros(matrizPlano,Tx,xn,yn);
Lmuros=Lmuros.*Lw;
%4 Calculo para suelos 
Lsuelo=Kf*Lf; % NO se concidera suelos 
%5 Calculo potencia recibida en Rx multi-wall model MWM
Pmwm=Pt+Gt+Gr-(Lfs+Lc+Lmuros+Lsuelo); 
%Graficas
figure(1)
hold on
Pwalls=coreccionPuntos(Pmwm(p,:));
plot(x,Pwalls);
legend('Espacio libre','Multi-Wall Model');
%Grafica en el plano 2D 
Pmwm=matrizCor(Pmwm)
x=resolucion:resolucion:dmax;
x=[-flip(x) 0 x];
figure(6)
hold on
PRgraf=imresize(mat2gray(Pmwm),[1000 1000]);
RI = imref2d(size(PRgraf)); RI.XWorldLimits = [(-dmax+Tx(1)) dmax+Tx(1)]; RI.YWorldLimits = [-dmax+Tx(2) dmax+Tx(2)];
imshow(PRgraf,RI);
%image(PRgraf,'CDataMapping','scaled');
colormap(gca,'jet');

grid on
hold on
min(Pmwm);
end

function [x,y]=limites(Vlineas)
%Grafica y puntos extremos: 
k=0;
x=0;y=0;
for i=1:length(Vlineas)
    %Geometrias  
    p=Vlineas{i};
    figure(2)
    hold on
    %Se toman solo las coordenadas (x,y)
    plot([p(1) p(4)],[p(2) p(5)]);
    title('Plano');
    hold off
    %Bordes.
    %Se almacenan los puntos de las geometrias. 
    x(i+k)=p(1);
    x(i+k+1)=p(4);
    y(i+k)=p(2);
    y(i+k+1)=p(5);
    k=k+1;
end   
end

%%
%{
 Funcion ubicarPuntos:
 Argumentos de entrada: matrizPlano, matriz de zeros para almacenar los
 puntos de las paredes del plano 
 x,y, vectores de distancia (Saltos en eje x y y) 
px,py (Coordenadas de las rectas)
%}
function matrizPlano=ubicarPuntos(matrizPlano,x,y,px,py)
    for n=1:2:length(px)
        %rectas verticales, coordenada 'x' no varia 
        if(px(n)==px(n+1))
            i=buscar(x,px(n));
            a=buscar(y,py(n));
            b=buscar(y,py(n+1));
            [a,b]=ordenar(a,b);
            %cordenadas (i,a) y (i,b). extermos del segmento
            %Ubicar matriz 
            for k=a:b
                matrizPlano(k,i)=1;
            end
        else 
            %rectas horizontales, coordenada 'y' no varia
            if(py(n)==py(n+1))
                i=buscar(y,py(n));
                a=buscar(x,px(n));
                b=buscar(x,px(n+1));
                [a,b]=ordenar(a,b);
                %cordenadas (i,a) y (i,b). extermos del segmento
                %Ubicar matriz 
                for k=a:b
                    matrizPlano(i,k)=1;
                end
            else 
            %rectas oblicuas, pendiente           
                disp('Recta oblicua pendiente');
            end
            
        end 
    end
end
%%
%{
 Funcion ordenar: ordena dos numeros de mayor a menor 
 Argumentos de entrada:a,b
salida: iinicio ifin 
%}
function [inicio, fin]=ordenar (a,b)
    inicio=min(a,b);
    fin=max(a,b);
end
%%
%{
 Funcion buscar: ubicar coordenada de la componente 
 Argumentos de entrada:x,p
salida: iinicio ifin 
%}
function [c]=buscar(x,p)
i=1;c=1e8;
    while (i<=length(x))
        if x(i)==p
            c=i;
            i=length(x)+1;
        else
            i=i+1;
        end
    end
end

%{
Funcion DeteccionMuros: Calculo La perdida en los muros segun la ubicacion
de la antena
Entrada: matrizPlano m (matriz de distancias, Tx (ubucacion de transmisor)
%} 
function matriz=DeteccionMuros(matrizPlano, Tx,xn,yn)
    %Ecuacion de la recta 
    yaux=flip(yn);
    matriz =matrizPlano;
    xt=Tx(1);yt=Tx(2);
    paso=xn(2)-xn(1);
    %Sector 1 .. ver en papel ycont=0 varia x
    for i=1:length(xn) 
        %Ecuacion de la recta  
        if(xt==max(xt,xn(i)))
            x1=xn(i);y1=yn(1);
            x2=xt;y2=yt;
        else
            x2=xn(i);y2=yn(1);
            x1=xt;y1=yt;
        end
        m=round((y2-y1)/(x2-x1),2) ;  
        b=-m*x1+y1;
        x=x1:paso:x2;
        fx=round(m*x+b,2);
        fx=aproximarVector(fx,paso);
        if (abs(m)>10)
           fx=y1:paso:y2;
           x=ones(length(fx)).*x1;
        end
        if (m<-10)
           a=find(xn==x1); 
           yv=matrizPlano(:,a);
           b=find(yn==y1);
           yv=flip(yv);
           yv=yv(1:b);
           b=find(yv==1);
           fx=y2:paso:yn(b);
           x=ones(length(fx)).*x2;
        end
        %Cubrir puntos 
        contmuros=0;
        if(m>=0)
            for k=1:length(x)
                a=find(xn==x(length(x)+1-k));
                b=find(yaux==fx(length(x)+1-k));
                contmuros=cruce(matrizPlano,b,a,contmuros);
                matriz(b,a)=contmuros;
            end
        else 
            for k=1:length(x)
                a=find(xn==x(k));
                b=find(yaux==fx(k));
                contmuros=cruce(matrizPlano,b,a,contmuros);
                matriz(b,a)=contmuros;
            end
        end
        
    end
     %Sector 2 ..ycont=yfinal  varia  x
    %%
    for i=1:length(xn) 
        %Ecuacion de la recta  
        if(xt==max(xt,xn(i)))
            x1=xn(i);y1=yn(length(yn));
            x2=xt;y2=yt;
        else
            x2=xn(i);y2=yn(length(yn));
            x1=xt;y1=yt;
        end
        m=round((y2-y1)/(x2-x1),2) ;  
        b=-m*x1+y1;
        x=x1:paso:x2;
        fx=round(m*x+b,2);
        fx=aproximarVector(fx,paso);
        if (abs(m)>10)
           fx=y1:paso:y2;
           x=ones(length(fx)).*x2;
        end
        if (m<-10)
           a=find(xn==x1); 
           yv=matrizPlano(:,a);
           b=find(yn==y1);
           yv=flip(yv);
           yv=yv(b:length(yv));
           b=find(yv==1);
           fx=yn(b):paso:y2;
           x=ones(length(fx)).*x1;
        end
        %Cubrir puntos 
        contmuros=0;
        if(m<=0)
            for k=1:length(x)
                a=find(xn==x(length(x)+1-k));
                b=find(yaux==fx(length(x)+1-k));
                contmuros=cruce(matrizPlano,b,a,contmuros);
                matriz(b,a)=contmuros;
            end
        else 
            for k=1:length(x)
                a=find(xn==x(k));
                b=find(yaux==fx(k));
                contmuros=cruce(matrizPlano,b,a,contmuros);
                matriz(b,a)=contmuros;
            end
        end
        
    end
    %%
    %Sector 3 x=xinicial y contante 
    for i=1:length(yn) 
        %Ecuacion de la recta  
        x1=xn(1);y1=yn(i);
        x2=xt;y2=yt;
        m=round((y2-y1)/(x2-x1),2) ;  
        b=-m*x1+y1;
        x=x1:paso:x2;
        fx=round(m*x+b,2);
        fx=aproximarVector(fx,paso);
        %Cubrir puntos 
        contmuros=0;
            for k=1:length(x)
                a=find(xn==x(length(x)+1-k));
                b=find(yaux==fx(length(x)+1-k));
                contmuros=cruce(matrizPlano,b,a,contmuros);
                matriz(b,a)=contmuros;
            end
    end
   
    %%
      %Sector 4 x1 = xt y conts
    for i=1:length(yn) 
        %Ecuacion de la recta  
        x2=xn(length(xn));y2=yn(i);
        x1=xt;y1=yt;
        m=round((y2-y1)/(x2-x1),2) ;  
        b=-m*x1+y1;
        x=x1:paso:x2;
        fx=round(m*x+b,2);
        fx=aproximarVector(fx,paso);
        %Cubrir puntos 
        contmuros=0;
            for k=1:length(x)
                a=find(xn==x(k));
                b=find(yaux==fx(k));
                contmuros=cruce(matrizPlano,b,a,contmuros);
                matriz(b,a)=contmuros;
            end
    end  
    figure(5)
    imshow(matriz);
    title('matriz');
end 
function c=cruce(m,b,a,mu)
    c=mu;
    if b>1 
        
    end 
    if(m(b,a)==1) 
        c=c+1;
    end
end
function v=aproximarVector(vector,paso)
    for i=1:length(vector)
        a=round(vector(i),1);
        a=vector(i)-a;
        dif=a-paso;
        if(abs(dif)<0.03)
            v(i)=vector(i)-dif;
        else     
            v(i)=round(vector(i),1);
        end
    end 
end
function Pmwm=coreccionPuntos(Pmwm)
vector=Pmwm;
pmaximo=find(vector==max(vector));
l=length(vector);
%Izquierda
for j=2:pmaximo-1
    if(vector(pmaximo-j)>vector(pmaximo -j+1))
        alpha=vector(pmaximo-j)-vector(pmaximo-j+1);
        vector(pmaximo-j)=vector(pmaximo -j)-alpha;
    end
end
%derecho 
for j=2:l-pmaximo-1
    if(vector(pmaximo+j)>vector(pmaximo+j-1))
        alpha=vector(pmaximo+j)-vector(pmaximo+j-1);
        vector(pmaximo+j)=vector(pmaximo+j)-alpha;
    end
end
Pmwm=vector;
end

function h=matrizCor(h)
[n r]=size(h);
for i=1:n
    h(i,:)=coreccionPuntos(h(i,:));
end
end

