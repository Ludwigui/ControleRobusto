clc
clearvars
yalmip('clear');
format long e
%%
% Inicializando variáveis

ps = 60*10^5;
pt = 0;
betae = 14*10^8;%Módulo de compressibilidade
Ba = 500; %Coeficiente de atrito viscoso
Ks1 = 50; %Ganho do transdutor
L = 0.2; %metros
Me = 100; %Massa do sistema
Kv = 1.43;
KqUo = 2*10^(-5);
KpUo = 4.7*10^(7);
KcUo = KqUo/KpUo;
Kx = 5982;
wnv = 377;
qsiv = 0.7;
areaembolo = 4.91*10^(-4);%mm^2
areahaste = 2.54*10^(-4);%mm^2
Aa = areaembolo-areahaste; 
Vt = Aa*L;

%% Função de transferência

s = tf('s');
G = KqUo/(s^3*(Vt*Me/(4*betae*Aa)) + s^2*(Vt*Ba/(4*betae*Aa)+KcUo*Me/Aa) + s*(Aa+Ba*KcUo/Aa+Vt*Kx/(4*betae*Aa)) + Kx*KcUo/Aa);
% t95% de malha aberta igual a 64.3 segundos.
%ganho de malha aberta:
gain = dcgain(G);
figure
step(G)
title('Sistema em malha aberta')

spacestates = ss(G);
%% Definindo o sistema
% A é uma matriz 3x3. Separarei em elementos para limpar o código.

a11=0; a12=1; a13=0;
a21=0; a22=0; a23=1;
a31=-Kx*KcUo/Aa*(4*betae*Aa)/(Vt*Me); a32= -(Aa + (Ba*KcUo/Aa) + Vt*Kx/(4*betae*Aa))*(4*betae*Aa)/(Vt*Me); a33 = -((Vt*Ba)/(4*betae*Aa) + Me*KcUo/Aa)*(4*betae*Aa)/(Vt*Me);

PS = [0.02; -0.02]; %Supondo que exista um erro de estimação de +/-2% na constante da mola.
A0 = [a11 a12 a13;
     a21 a22 a23;
     a31 a32 a33];
A1 = [0 0 0 0;
      0 0 0 0;
      -Kx*KcUo*4*betae/(Vt*Me) -Kx/Me 0 0;
      0 0 0 0];

  

B = [0;0;KqUo*(4*betae*Aa)/(Vt*Me)];
Bw = [0;0;-KcUo/Aa*(4*betae*Aa)/(Vt*Me)];
C = [1 0 0];
D = 0;

%matrizes aumentadas
AA = [A0 zeros(3,1);
    -C 0];
Bua = [B; 0];
Bwa = [Bw; 0];
Ca = [C 0];

Ac = 1;
Bc = 1;
Cc = 1;

%% Definição de malha fechada:
  
%%
nx = size(AA,1);
%nv = size(PS,1);
nr = size(Bua,2);
nq = size(Ca,1);
np = size(Bwa,2);


[nl,nv] = size(PS);

% Decision variables
P = sdpvar(nx,nx);
W = sdpvar(nr,nx,'full');
M = sdpvar(nq,nq);
g2 = sdpvar(1,1);
options = sdpsettings('verbose',0,'solver','sdpt3');
%% polos com parte real menor que 3
 L1=6;
 M1=1;
% %% dentro de uma circunferência de raio 4
L2=[-4 0;0 -4];
M2=[0 1;0 0];
% %% abaixo de determinado ângulo oriundo da origem do plano complexo
% M3 = 0.5*sqrt(2)*[1 1; -1 1];
% nm = size(M3,1);
% L3 = zeros(nm);
%% LMIs

% LMIs definition
LMIs = [g2 - trace(M) >= 0] + [[M (Ca*P+D*W); (Ca*P+D*W)' P] >= 0];
for i=1:nv
    Ad = AA+PS(i)*A1;
    LMIs = LMIs +  [(Ad*P + Bua*W)+(Ad*P+Bua*W)' + Bwa*Bwa' <= 0] + [kron(L1,P) + kron(M1,(Ad*P+Bua*W)) + kron(M1',(Ad*P+Bua*W)') <= 0] + [kron(L2,P) + kron(M2,(Ad*P+Bua*W)) + kron(M2',(Ad*P+Bua*W)') <= 0];
end
%%
% Calling the solver
result = optimize(LMIs,g2,options);
%--------------------------------------------------------------------------
% Recovering the results
P = double(P);
W = double(W);
K = W*inv(P);
g2 = double(g2);
g2 = sqrt(g2);

%--------------------------------------------------------------------------
% Result validation
test_LMI = check(LMIs);
if test_LMI >= 0
    disp('The solution is reliable !!!');
    disp(result);
    disp(test_LMI);
    disp('The uncertain system has feasible solution:');
    display(P);
    display(K);
    
else
    disp('The solution is not reliable :(');
    msg = sprintf('%s codigo nro: %d',result.info,result.problem);
    disp(msg);
    disp(result);
    disp(test_LMI);
end

