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
wnv = 377;
qsiv = 0.7;
areaembolo = 4.91*10^(-4);%mm^2
areahaste = 2.54*10^(-4);%mm^2
Aa = areaembolo-areahaste; 
Vt = Aa*L;
%% Definindo o sistema
% A é uma matriz 3x3. Separarei em elementos para limpar o código.

a11=0; a12=1; a13=0;
a21=0; a22=0; a23=1;
a31=0; a32= -(Aa+(Ba/Aa)*KcUo)*4*betae*Aa/(Vt*Me); a33 = -(Vt/(4*betae)*(Ba/Aa)+KcUo*(Me/Aa))*4*betae*Aa/(Vt*Me);

A = [a11 a12 a13;
     a21 a22 a23;
     a31 a32 a33];
B = [0;0;KqUo];
%% Hurwitz stability

nx = size(A,1);
% Variável de decisão
P = sdpvar(nx,nx);
options = sdpsettings('verbose',0,'solver','sdpt3');
% Definição de LMIs

LMIs = [];
LMIs = LMIs + [P>=1e-3] + [A'*P + P*A <=0];

result = optimize(LMIs,[],options);
P = double(P);

% Result validation
test_LMI = check(LMIs);
if test_LMI >= 0
    disp('The solution is reliable !!!');
    disp(result);
    disp(test_LMI);
    disp('A feasible solution is P =');
    display(P);
    
else
    disp('The solution is not reliable :');
    msg = sprintf('%s codigo nro: %d',result.info,result.problem);
    disp(msg);
    disp(result);
    disp(test_LMI);
end
%% sistema em malha aberta

[t,x] = ode45(@sysopenloop,[0 2],[1 0 0]);
plot(t,x)
grid
legend('x_1','x_2','x_3')
title('Hurwitz Estável?')

%% 