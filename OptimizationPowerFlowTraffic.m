%        Trabalho realizado por Alexandre Amorim N�1161497

%----------------------------------AlgoritmoResolu��oRedes---------------------------------------

%% Introdu��o de dados

clc
clear all
     
% Dados de entrada do programa

% SBase=input('Insira a pot�ncia aparente de base - ');
% linhas=input ('Insira a matriz linhas - ');
% cargas=input('Insira a matriz cargas - ');
% prod=input('Insira a matriz produ��o - ');
% tipo=input('Insira a matriz tipo - ');

    % Rede com 3 barramentos

linhas=[1 2 0.01 0.1 0.001; 1 3 0.01 0.1 0.001; 2 3 0.02 0.2 0.002];
cargas=[1 0.7 0.2; 2 1 0.3; 3 0.5 0.1];
prod=[2 0.9 1.01; 3 1 1.02];
tipo=['Q','V','R']; %Q = Barramento PQ; P = Barramento PV; R = Barramento de Refer�ncia dos Argumentos 

    % Rede com 4 Barramentos
    
%linhas=[1 2 0.01008 0.0504 0.1025; 1 3 0.00744 0.0372 0.0775; 2 4 0.00744 0.0372 0.0775; 3 4 0.01272 0.0636 0.06375];
%cargas=[1 0.5 0.3099; 2 1.7 1.0535; 3 2.6 1.2394; 4 0.8 4.958];
%prod=[1 1 1; 4 3.18 1.02];
%tipo=['R','Q','Q','V'];

    % Rede com 14 barramentos 

%linhas=[1 2 0.01938 0.05917 0.0528;1 5 0.05403 0.22304 0.0492;2 3 0.04699 0.19797 0.0438;2 4 0.05811 0.17632 0.0340;2 5 0.05695 0.17388 0.0346;3 4 0.06701 0.17103 0.0128;4 5 0.01335 0.04211 0;4 7 0.0 0.20912 0;4 9 0.0 0.55618 0;5 6 0.0 0.25202 0;6 11 0.09498 0.19890 0;6 12 0.12291 0.25581 0;6 13 0.06615 0.13027 0;7 8 0.0 0.17615 0;7 9 0.0 0.11001 0;9 10 0.03181 0.08450 0;9 14 0.12711 0.27038 0;10 11 0.08205 0.19207 0;12 13 0.22092 0.19988 0;13 14 0.17093 0.34802 0];
%cargas=[1 0 0;2 0.217 0.127;3 0.942 0.19;4 0.478 -0.039;5 0.076 0.016;6 0.112 0.075;7 0 0;8 0 0;9 0.295 0.166;10 0.090 0.058;11 0.035 0.018 ;12 0.061 0.016;13 0.135 0.058;14 0.149 0.050];
%prod=[1 2.324 1.060;2 0.4 1.045;3 0 1.010;6 0 1.070;8 0 1.090]; 
%tipo=['R','V','V','Q','Q','V','Q','V','Q','Q','Q','Q','Q','Q'];

% Menu para escolha do m�todo pretendido
% Mostra os diferentes m�todos que podem ser selecionados

     disp('Escolha o m�todo pretendido')
     disp(' ')
     disp('1 -> Para Modelo DC')
     disp('2 -> Para M�todo de Gauss-Seidel')
     disp('3 -> Para M�todo de Newton-Raphson')
     disp('4 -> Para M�todo R�pido Desacoplado')
     disp(' ')
     metodo=input('Introduza o m�todo pretendido: '); % Escolha do m�todo a executar
     
     %Eliminar poss�veis erros de escolha do utilizador
     
     if metodo<1 || metodo>4 
        clc
        disp('O n�mero que introduziu � inv�lido.') % Erro em caso de uma entrada inv�lida
        disp('Por favor, tente novamente.')
        disp(' ')
        disp('1 -> Para Modelo DC')
        disp('2 -> Para M�todo de Gauss-Seidel')
        disp('3 -> Para M�todo de Newton-Raphson')
        disp('4 -> Para M�todo de R�pido Desacoplado')
        disp(' ')
        metodo=input('Introduza o m�todo pretendido: ');
     else
        clc
     end
     
     
%% C�lculos iniciais 

Origem=linhas(:,1);           %N� de Origem
Destino=linhas(:,2);          %N� de Destino
R=linhas(:,3);                %Resist�ncias das linhas
X=linhas(:,4);                %Reat�ncias das linhas
Ysh=(linhas(:,5))*1j;         %Admit�ncias shunt das linhas
Barramento=1:length(tipo);    %Id dos barramentos
Pc=cargas(:,2);               %Pot�ncia ativa consumida pelas cargas         
Pg=zeros(length(tipo),1);     %Pot�ncia ativa gerada no barramento'
P1=prod(:,1);                 %Pot�ncia ativa gerada no barramento
Z=R+j*X;                      %Imped�ncias das linhas
Qc=cargas(:,3);               %Pot�ncia Reativa consumidas pelas cargas
Qg=zeros(length(tipo),1);     %Pot�ncia Reativa gerada no barramento'
tQc=Qc';                      %Pot�ncia Reativa consumida no barramento'
    
%% Sele��o do m�todo de c�lculo

switch metodo 
    
%% Modelo DC 

    case 1
%if R < X  (na verdade � R << X, logo R ~ 0)      
for i=1:length(tipo)
    for ii=1:length(prod(:,1))
    if Barramento(i)==prod(ii,1)
        Pg(i)=prod(ii,2);
    end
    end
end

P=Pg-Pc;                                           %Pot�ncia injetada no n�

tic;                                           %In�cio da contagem do tempo

%C�lculo da matriz das suscept�ncias
B=zeros(length(tipo));
for i=1:length(Origem)
       B((linhas(i,1)),(linhas(i,2)))=(-1/(linhas(i,4)));     %Bik = -1/Xik
       B((linhas(i,2)),(linhas(i,1)))=(-1/(linhas(i,4)));     %Bki = -1/Xki    
end 

for i=1:length(tipo)
    h=0;
    for ii=1:length(tipo) 
    h=h+abs(B((i),(ii)));                       %Bii= S(-1/Xik).-1
    end
    B((i),(i))=h;
end
for i=1:length(tipo)
    if tipo(i)=='R'        %Se for o BarRef eliminamos a sua linha e coluna      
        B(i,:)=[];
        B(:,i)=[];
        P(i,:)=[];
        BarramentoRef=i;
    end
end

t=(B^-1)*P;                           %Vetor coluna dos �ngulos das tens�es

for i=1:length(tipo)
    if i<BarramentoRef     %Se o Id do barramento for inferior ao do BarRef
        teta(i,1)=t(i);
    elseif i==BarramentoRef   %Se o Id do barramento for igual ao do BarRef
        teta(i,1)=0;
    else
        teta(i,1)=t(i-1);  %Se o Id do barramento for superior ao do BarRef 
    end
end

for i=1:length(Origem)
    FluxoPot(i,1)=(teta(Origem(i))-teta(Destino(i)))/(X(i));   
end

Pij(1,:)=Origem;                                              %N� do Origem
Pij(2,:)=Destino;                                            %N� de Destino
Pij(3,:)=FluxoPot;

Time=toc;                                         %Fim do tempo de contagem

disp('-------------------------Modelo DC----------------------------')
disp(' ')
display(['Tempo ',num2str(Time),' segundos']);
disp(' ')
for i=1:length(tipo)  
    disp(['teta ',num2str(i),' = ',num2str(teta(i),'%.4f')]);
end
disp(' ')
fprintf('Pij =\n\n')
fprintf('    %d    %d    %.4f\n',Pij)

%% M�todo Gauss-Seidel

 case 2
        
for i=1:length(tipo)
    for ii=1:length(prod(:,1))
    if Barramento(i)==prod(ii,1)
        Pg(i)=prod(ii,2);
    end
    end
end
%--------------------------------------------------------------------------
% C�lculo da matriz y
Y=zeros(length(tipo));
for i=1:length(Origem)
    Y(Origem(i),Destino(i))=-(1/Z(i));                                     %Yik = -1/Zik
    Y(Destino(i),Origem(i))=-(1/Z(i));                                     %Yki = -1/Zki
    Y(Origem(i),Origem(i))=Y(Origem(i),Origem(i))+ 1/Z(i)+Ysh(i)/2;%Yii = S(-1/Zik)+Ysh(ik)/2
    Y(Destino(i),Destino(i))=Y(Destino(i),Destino(i))+ 1/Z(i)+Ysh(i)/2;%Ykk = S(-1/Zki)+Ysh(ki)/2
end

%Solu��o Inicial
itermax=200;
erromax=0.01; %Erro m�ximo definido para os desvios de pot�ncia
U=zeros(length(tipo),1);
U(P1(:))=prod(:,3);
q=zeros(length(tipo),1);
p=zeros(length(tipo),1);

    for i=1:length(tipo)
        if tipo(i)=='Q'
            %U(i) = 1; Calcular P e Q para barramento PQ
            U(i)=1;
            p(i,1)=Pg(i)-Pc(i);
            q(i,1)=Qg(i)-Qc(i);
        elseif tipo(i)=='V'
            %Calcular P para barramento PV
            p(i,1)=Pg(i)-Pc(i);
        end
    end
Uini=U;
   
disp('-----------------M�todo Gauss Seidel------------------');
tic;                                                    %In�cio da contagem
flag=1;
niter=0;                                               %N�mero de Itera��es
while flag == 1
    Ical=zeros(length(tipo),1);
    for i=1:length(tipo)
        if tipo(i)=='V'                                      %Barramento PV
            for ii=1:length(tipo)
                %C�lculo da Matriz das correntes
                Ical(i)=Ical(i)+(Y(i,ii)*U(ii));         %Ical(i)=S(Yik.Uk)
            end
            Qcal(i)=imag(U(i)*conj(Ical(i)));         %Scal=Pcal+jQcal=U.I*
        end
    end   

    for i=1:length(tipo)
        if tipo(i)=='Q'                                      %Barramento PQ
            c=0;
            for ii=1:length(tipo)
                if i~=ii                             %Eliminar Yii (Yik.Ui)
                    c=c+(Y(i,ii)*U(ii));
                end
            end
            %Novas tens�es
            %C�lculo da tens�o para PQ                               %1  [P(i)-jQ(i)
            U(i)=(1/Y(i,i))*[(p(i)-1j*(q(i)))/(conj(U(i)))-c]; %U(i)=---.[---------  - S(Yik.Ui)]
        elseif tipo(i)=='V'                                         %Yii [  U(i)*               ]
            c=0;
            for ii=1:length(tipo)
                if i~=ii                             %Eliminar Yii S(Yik.Ui)
                    c=c+(Y(i,ii)*U(ii));
                end
            end
            %C�lculo e corre��o da tens�o para PV
            U(i)=(1/Y(i,i))*[(p(i)-1j*(Qcal(i)))/(conj(U(i)))-c];       %Uini(i)
            U(i)=U(i)*Uini(i)/abs(U(i));                          %U(i)=--------   
        end                                                             %|U(i)|
    end
% C�lculo de S para Obter Pcal e Qcal
    Scal=zeros(length(tipo),1);

    for i=1:length(tipo)      
        if tipo(i)~='R'                       %N�o se calcula para o BarRef
           Ical(i)=(Y(i,:)*U);                              %Ical=S(Yik.Uk)
        end
        Scal(i)=U(i)*conj(Ical(i));                     %Scal(i)=U(i).I(i)*
    end   

    Pcal=zeros(length(tipo),1);
    Qcal=zeros(length(tipo),1);
    for i=1:(length(tipo));
        Pcal(i)=real(Scal(i));
        Qcal(i)=imag(Scal(i));
    end
%Verifica��o de converg�ncia
    for i=1:length(tipo)
        if tipo(i)=='V'
            dp(i)=abs(p(i)-Pcal(i)); %dp = deltaP = Desvio das Pot�ncias Ativas
        elseif tipo(i)=='Q'
            dp(i)=abs(p(i)-Pcal(i)); %dp = deltaP = Desvio das Pot�ncias Ativas
            dq(i)=abs(q(i)-Qcal(i)); %dq = deltaQ = Desvio das Pot�ncias Reativas
        end
    end
    %Caso n�o convergiu => flag continua a 1; ciclo continua; niter += 1
    %Itera��o/Incrementar a vari�vel de n�mero de itera��es
    niter=niter+1;
    if max(dp)<0.01 & max(dq)<0.01
        %Se convergiu Flag a 0 e acaba
        flag=0;
        disp(' ')
        disp(['GS -> ',num2str(niter),' itera��es']);
        disp(' ')
        Time=toc;                                      %Fim da contagem do tempo
display(['Tempo ',num2str(Time),' segundos']);
disp(' ')
    end
    if niter>itermax
        flag=0;
        disp('N�o Converge');
    end
end
%Se tipo do barramento ~= Q => Calcular Ical e Scal
ical=zeros(length(prod(:,1)));
    for i=1:length(tipo)
        if tipo(i)~='Q'
            for ii=1:length(tipo)
                ical(i)=ical(i)+(Y(i,ii)*U(ii));
                Qcal(i)=imag(U(i)*conj(ical(i)));
                Qg(i)=Qcal(i)+tQc(i);                 %Qg(i)=Qcal(i)+Qc(i)'
            end
        end
    end
tQg=Qg';
for i=1:length(tipo)  
    disp(['U',num2str(i),' = ',num2str(U(i),'%.4f')]);
end
disp(' ')
%Manter m�dulo e alterar o argumento
for i=1:length(tipo)   
    Res = abs(U(i));
    teta = atan2(imag(U(i)),real(U(i)));
    disp(['U',num2str(i),' = ',sprintf('%.4f%c%.4f',Res,' ',teta)]);
end

for i=1:length(Origem)
    Z(Origem(i),Destino(i))=linhas(i,3)+j*linhas(i,4);
end
disp(' ');
%disp('Sik(Fluxo de Potencia de I para k)=');
for i=1:length(Origem)
    %Sik=(U(i).(U(i)-U(k))*)/Zik
    S(Origem(i),Destino(i))=U(Origem(i))*conj((U(Origem(i))-U(Destino(i)))/Z(Origem(i),Destino(i)));
    %Ski=(U(k).(U(k)-U(i))*)/Zki
    S(Destino(i),Origem(i))=U(Destino(i))*conj((U(Destino(i))-U(Origem(i)))/Z(Origem(i),Destino(i)));
    %Perdasik=Sik+Ski
    perdas(Origem(i),Destino(i))=S(Origem(i),Destino(i))+S(Destino(i),Origem(i));
end
S
disp('Perdas (Linha ik)');
perdas
tQg=tQg(prod(:,1));
Qgerado(1,:)=P1;
Qgerado(2,:)=tQg;                                       %Qgerado calculado'
disp('Qgerado');
fprintf('  %d    %.5f\n',Qgerado)

%% M�todo Newton-Raphson

    case 3
        
for i=1:length(tipo)
    for ii=1:length(prod(:,1))
    if bar(i)==prod(ii,1);
        Pg(i)=prod(ii,2);
    end
    end
end
%-----------------------------------------------------------------------

%Matriz y
Y=zeros(length(tipo));
for i=1:length(Origem)
    Y(Origem(i),Destino(i))=-(1/Z(i));
    Y(Destino(i),Origem(i))=-(1/Z(i));
    Y(Origem(i),Origem(i))=Y(Origem(i),Origem(i))+ 1/Z(i)+Ysh(i)/2;
    Y(Destino(i),Destino(i))=Y(Destino(i),Destino(i))+ 1/Z(i)+Ysh(i)/2;
end

%Solu��o inicial
itermax=100;
erromax=0.01;
U=zeros(length(tipo),1);
U(P1(:))=prod(:,3);
Q=zeros(length(tipo),1);
P=zeros(length(tipo),1);

%Calculo das pot�ncias injetadas nos barramentos PV e PQ
    for i=1:length(tipo)
        if tipo(i)=='Q'
            U(i)=1;
            P(i,1)=Pg(i)-Pc(i);
            Q(i,1)=Qg(i)-Qc(i);
        elseif tipo(i)=='V'
            P(i,1)=Pg(i)-Pc(i);
        end
    end
Uini=U;
tetaini=zeros(length(Uini),1);
   
disp('--------------------M�todo Newton-Raphson--------------------');
disp(' ' )
tic;
flag=1;
niter=0;
A=find(tipo==('V'));
B=find(tipo==('Q'));
C(:,1)=[A,B];
C=sort(C);                                          %ordena ascendentemente
%C�lculo das potencias para barramentos PV
    Ical=zeros(length(tipo),1);   
    for i=1:length(tipo)
        if tipo(i)=='V'
            for ii=1:length(tipo)
                Ical(i)=Ical(i)+(Y(i,ii)*Uini(ii));
            end
            Scal(i)=(Uini(i)*conj(Ical(i)));
            Pcal(i,1)=real(Scal(i));
            Qcal(i,1)=imag(Scal(i));
        end
    end  
%C�lculo das potencias para barramentos PQ
    for i=1:length(tipo)
        if tipo(i)=='Q'
            for ii=1:length(tipo)
                Ical(i)=Ical(i)+(Y(i,ii)*Uini(ii));
             end
            Scal(i)=(Uini(i)*conj(Ical(i)));
            Pcal(i,1)=real(Scal(i));
            Qcal(i,1)=imag(Scal(i));
         end
    end 
    
v=abs(Uini);

while flag == 1
    dp1=zeros(length(tipo),1);
    dq1=zeros(length(tipo),1);
        for i=1:length(A)
            dp1(A(i),1)=P(A(i))-Pcal(A(i));
        end
        for i=1:length(B)
            dp1(B(i),1)=P(B(i))-Pcal(B(i));
            dq1(B(i),1)=Q(B(i))-Qcal(B(i));
        end
        
clear dd1
dd1=[dp1(C);dq1(B)];
%Preenchimento da matriz Jacobiana
H=zeros(length(tipo)-1);
M=zeros(length(tipo)-length(prod(:,1)),length(tipo)-1);
N=zeros(length(tipo)-1,length(tipo)-length(prod(:,1)));
L=zeros(length(tipo)-length(prod(:,1)),length(tipo)-length(prod(:,1)));

for i=1:length(H(:,1))
    for ii=1:length(H(1,:))
        if i==ii
            if C(i)~=A
            H(i,ii)=-Q(C(i),1)-imag(Y(C(i),C(i)))*v(C(i))^2;
            else
            H(i,ii)=-Qcal(C(i),1)-imag(Y(C(i),C(i)))*v(C(i))^2;
            end
        else
            H(i,ii)=v(C(i))*v(C(ii))*(real(Y(C(i),C(ii)))*sin(tetaini(C(i))-tetaini(C(ii)))-imag(Y(C(i),C(ii)))*cos(tetaini(C(i))-tetaini(C(ii))));
        end  
    end
end

for i=1:length(L(:,1))
    for ii=1:length(L(1,:))
            if i==ii
                if C(i)~=A
                  L(i,ii)=Q(C(i),1)-imag(Y(C(i),C(i)))*v(C(i))^2;
                else
                    L(i,ii)=Qcal(C(i),1)-imag(Y(C(i),C(i)))*v(C(i))^2;
                end
                else
                    L(i,ii)=v(C(i))*v(C(ii))*(real(Y(C(i),C(ii)))*sin(tetaini(C(i))-tetaini(C(ii)))-imag(Y(C(i),C(ii)))*cos(tetaini(C(i))-tetaini(C(ii))));
            end
    end
end

for i=1:length(N(:,1))
    for ii=1:length(N(1,:))
            if i==ii
                  N(i,ii)=P(C(i),1)+real(Y(C(i),C(ii)))*v(C(i))^2;
                else
              N(i,ii)=v(C(i))*v(C(ii))*(real(Y(C(i),C(ii)))*cos(tetaini(C(i))-tetaini(C(ii)))+imag(Y(C(i),C(ii)))*sin(tetaini(C(i))-tetaini(C(ii))));
            end
    end
end

for i=1:length(M(:,1))
    for ii=1:length(M(1,:))
            if i==ii
                  M(i,ii)=P(C(i),1)-real(Y(C(i),C(ii)))*v(C(i))^2;
                else
              M(i,ii)=-v(C(i))*v(C(ii))*(real(Y(C(i),C(ii)))*cos(tetaini(C(i))-tetaini(C(ii)))+imag(Y(C(i),C(ii)))*sin(tetaini(C(i))-tetaini(C(ii))));
            end
    end
end

clear J1 J2 J t
J1=[H;M];
J2=[N;L];
J=[H,N;M,L];
t=inv(J)*dd1;
v=abs(Uini);

for i=1:length(B)
    v(B(i))=v(B(i))+t(length(C)+i);
end

for i=1:length(C)
    tetaini(C(i))=tetaini(C(i))+t(i);
end
[X,Z] = pol2cart(tetaini,v);
Uini=X+1i*Z;
    Ical=zeros(length(tipo),1); 
    clear Scal Pcal Qcal
    for i=1:length(tipo)
        if tipo(i)=='V'
            for ii=1:length(tipo)
                Ical(i)=Ical(i)+(Y(i,ii)*Uini(ii));
            end
            Scal(i)=(Uini(i)*conj(Ical(i)));
            Pcal(i,1)=real(Scal(i));
            Qcal(i,1)=imag(Scal(i));
        end
    end  
    for i=1:length(tipo)
        if tipo(i)=='Q'
            for ii=1:length(tipo)
                Ical(i)=Ical(i)+(Y(i,ii)*Uini(ii));
             end
            Scal(i)=(Uini(i)*conj(Ical(i)));
            Pcal(i,1)=real(Scal(i));
            Qcal(i,1)=imag(Scal(i));
         end
    end   
    clear dp dq
        for i=1:length(A)
            dp(A(i),1)=abs(P(A(i))-Pcal(A(i)));
        end
        for i=1:length(B)
            dp(B(i),1)=abs(P(B(i))-Pcal(B(i)));
            dq(B(i),1)=abs(Q(B(i))-Qcal(B(i)));
        end

    niter=niter+1;
    if max(dp)<0.01 & max(dq)<0.01
        flag=0;
        disp(['NR -> ',num2str(niter),' itera��es']);
    end
    if niter>itermax
        flag=0;
        disp('N�o Converge');
    end
end
Time=toc;
disp(' ')
display(['Tempo  ',num2str(Time),' segundos']);
disp(' ')
final=cell(length(v),3);
for i=1:length(v)
   final{i,1}=strcat('u',int2str(i));
   final{i,2}=num2str(v(i));
   final{i,3}=num2str(tetaini(i));
end
disp('Tens�es')
disp(v)
disp('Argumentos')
disp(tetaini)

for i=1:length(Origem)
    Z(Origem(i),Destino(i))=linhas(i,3)+1i*linhas(i,4);
end

perdas=zeros(length(Origem));
for i=1:length(Origem)
    S(Origem(i),Destino(i))=Uini(Origem(i))*conj((Uini(Origem(i))-Uini(Destino(i)))/Z(Origem(i),Destino(i)));
    S(Destino(i),Origem(i))=Uini(Destino(i))*conj((Uini(Destino(i))-Uini(Origem(i)))/Z(Origem(i),Destino(i)));
    perdas(Origem(i),Destino(i))=S(Origem(i),Destino(i))+S(Destino(i),Origem(i));
end
disp('Sik(Fluxo de Pot�ncia de i para k)');
disp(S);
disp('Perdas (Linha I k)');
disp(perdas)

ical=zeros(length(prod(:,1)));
    for i=1:length(tipo)
        if tipo(i)~='Q'
            for ii=1:length(tipo)
                ical(i)=ical(i)+(Y(i,ii)*Uini(ii));
                Qcal(i)=imag(Uini(i)*conj(ical(i)));
                Qg(i)=Qcal(i)+tQc(i);
            end
        end
    end
tQg=Qg';
tQg=tQg(prod(:,1));
Qgerado(1,:)=P1;
Qgerado(2,:)=tQg;

disp('Qgerado');
fprintf('  %d    %.5f\n',Qgerado)

%% M�todo R�pido Desacopolado

    case 4

for i=1:length(tipo)
    for ii=1:length(prod(:,1))
    if bar(i)==prod(ii,1);
        Pg(i)=prod(ii,2);
    end
    end
end
%-----------------------------------------------------------------------
%matriz y
Y=zeros(length(tipo));
for i=1:length(Origem)
    Y(Origem(i),Destino(i))=-(1/Z(i));
    Y(Destino(i),Origem(i))=-(1/Z(i));
    Y(Origem(i),Origem(i))=Y(Origem(i),Origem(i))+ 1/Z(i)+Ysh(i)/2;
    Y(Destino(i),Destino(i))=Y(Destino(i),Destino(i))+ 1/Z(i)+Ysh(i)/2;
end
%solu��o inicial
itermax=100;
erromax=0.01;
U=zeros(length(tipo),1);
U(P1(:))=prod(:,3);
q=zeros(length(tipo),1);
p=zeros(length(tipo),1);

    for i=1:length(tipo)
        if tipo(i)=='Q'
            U(i)=1;
            p(i,1)=Pg(i)-Pc(i);
            q(i,1)=Qg(i)-Qc(i);
        elseif tipo(i)=='V'
            p(i,1)=Pg(i)-Pc(i);
        end
    end
teta=zeros(length(U),1);
   
disp('M�todo R�pido Desacoplado');
tic;
flag=1;
niter=0;
A=find(tipo==('V'));
B=find(tipo==('Q'));
C(:,1)=[A,B];
C=sort(C);

    Ical=zeros(length(tipo),1);   
    for i=1:length(tipo)
        if tipo(i)=='V'
            for ii=1:length(tipo)
                Ical(i)=Ical(i)+(Y(i,ii)*U(ii));
            end
            Scal(i)=(U(i)*conj(Ical(i)));
            Pcal(i,1)=real(Scal(i));
            Qcal(i,1)=imag(Scal(i));
        end
    end  
    for i=1:length(tipo)
        if tipo(i)=='Q'
            for ii=1:length(tipo)
                Ical(i)=Ical(i)+(Y(i,ii)*U(ii));
             end
            Scal(i)=(U(i)*conj(Ical(i)));
            Pcal(i,1)=real(Scal(i));
            Qcal(i,1)=imag(Scal(i));
         end
    end 

while flag == 1
    clear dp1 dq1
    dp1=zeros(length(tipo),1);
    dq1=zeros(length(tipo),1);
        for i=1:length(A)
            dp1(A(i),1)=p(A(i))-Pcal(A(i));
        end
        for i=1:length(B)
            dp1(B(i),1)=p(B(i))-Pcal(B(i));
            dq1(B(i),1)=q(B(i))-Qcal(B(i));
        end
b1=zeros(length(tipo)-1);
b2=zeros(length(tipo)-length(prod(:,1)),length(tipo)-length(prod(:,1)));
for i=1:length(b1(:,1))
    for ii=1:length(b1(1,:))
        b1(i,ii)=-imag(Y(i,ii));
    end
end
for i=1:length(b2(:,1))
    for ii=1:length(b2(1,:))
    b2(i,ii)=-imag(Y(i,ii));
    end
end
dtheta=inv(b1)*dp1(C);
dv=inv(b2)*dq1(B);
v=abs(U);
for i=1:length(B)
    v(B(i))=v(B(i))+dv(i);
end

for i=1:length(C)
    teta(C(i))=teta(C(i))+dtheta(i);
end
[X,Z] = pol2cart(teta,v);
clear U
U=X+1i*Z;
    Ical=zeros(length(tipo),1); 
    clear Scal Pcal Qcal
    for i=1:length(tipo)
        if tipo(i)=='V'
            for ii=1:length(tipo)
                Ical(i)=Ical(i)+(Y(i,ii)*U(ii));
            end
            Scal(i)=(U(i)*conj(Ical(i)));
            Pcal(i,1)=real(Scal(i));
            Qcal(i,1)=imag(Scal(i));
        end
    end  
    for i=1:length(tipo)
        if tipo(i)=='Q'
            for ii=1:length(tipo)
                Ical(i)=Ical(i)+(Y(i,ii)*U(ii));
             end
            Scal(i)=(U(i)*conj(Ical(i)));
            Pcal(i,1)=real(Scal(i));
            Qcal(i,1)=imag(Scal(i));
         end
    end   
    clear dp dq
        for i=1:length(A)
            dp(A(i),1)=abs(p(A(i))-Pcal(A(i)));
        end
        for i=1:length(B)
            dp(B(i),1)=abs(p(B(i))-Pcal(B(i)));
            dq(B(i),1)=abs(q(B(i))-Qcal(B(i)));
        end

    niter=niter+1;
    if max(dp)<0.01 & max(dq)<0.01
        flag=0;
        disp(['FDL -> ',num2str(niter),' itera��es']);
    end
    if niter>itermax
        flag=0;
        disp('N�o Converge');
    end
end
Time=toc;
display(['Tempo  ',num2str(Time),' segundos']);
final=cell(length(v),3);
for i=1:length(v)
   final{i,1}=strcat('u',int2str(i));
   final{i,2}=num2str(v(i));
   final{i,3}=num2str(teta(i));
end
final
for i=1:length(Origem)
    Z(Origem(i),Destino(i))=linhas(i,3)+1i*linhas(i,4);
end
disp('Sik(Fluxo de Potencia de I para k)=');
for i=1:length(Origem)
    S(Origem(i),Destino(i))=U(Origem(i))*conj((U(Origem(i))-U(Destino(i)))/Z(Origem(i),Destino(i)));
    S(Destino(i),Origem(i))=U(Destino(i))*conj((U(Destino(i))-U(Origem(i)))/Z(Origem(i),Destino(i)));
    perdas(Origem(i),Destino(i))=S(Origem(i),Destino(i))+S(Destino(i),Origem(i));
end
S
disp('Perdas (Linha I k)');
perdas

ical=zeros(length(prod(:,1)));
    for i=1:length(tipo)
        if tipo(i)~='Q'
            for ii=1:length(tipo)
                ical(i)=ical(i)+(Y(i,ii)*U(ii));
                Qcal(i)=imag(U(i)*conj(ical(i)));
                Qg(i)=Qcal(i)+tQc(i);
            end
        end
    end
tQg=Qg';
tQg=tQg(prod(:,1));
Qgerado(1,:)=P1;
Qgerado(2,:)=tQg;
disp('Qgerado');
fprintf('  %d    %.5f\n',Qgerado)
end

