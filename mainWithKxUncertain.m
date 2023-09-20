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

%% Definindo o sistema
% A é uma matriz 3x3. Separarei em elementos para limpar o código.

a11=0; a12=1; a13=0;
a21=0; a22=0; a23=1;
a31=-Kx*KcUo/Aa*(4*betae*Aa)/(Vt*Me); a32= -(Aa + (Ba*KcUo/Aa) + Vt*Kx/(4*betae*Aa))*(4*betae*Aa)/(Vt*Me); a33 = -((Vt*Ba)/(4*betae*Aa) + Me*KcUo/Aa)*(4*betae*Aa)/(Vt*Me);

PS = [0.02; -0.02];

A0 = [a11 a12 a13;
     a21 a22 a23;
     a31 a32 a33];

A1 = [0 0 0;
      0 0 0;
      -Kx*KcUo*4*betae/(Vt*Me) -Kx/Me 0];

%%
nx = size(A0,1);
nv = size(PS,1);

% Decision variables
P = sdpvar(nx,nx);
options = sdpsettings('verbose',0,'solver','sdpt3');

B = [0;0;KqUo*(4*betae*Aa)/(Vt*Me)];

C = [1;0;0];

% LMIs definition
LMIs = [P >= 0];
for i=1:nv
    Ad = A0+PS(i)*A1;
    LMIs = LMIs + [Ad'*P+P*Ad <= 0];
end

%%
% Calling the solver
result = optimize(LMIs,[],options);
%--------------------------------------------------------------------------
% Recovering the results
P = double(P);
%--------------------------------------------------------------------------
% Result validation
test_LMI = check(LMIs);
if test_LMI >= 0
    disp('The solution is reliable !!!');
    disp(result);
    disp(test_LMI);
    disp('The uncertain system is quadratic stable with');
    display(P);
    
else
    disp('The solution is not reliable :(');
    msg = sprintf('%s codigo nro: %d',result.info,result.problem);
    disp(msg);
    disp(result);
    disp(test_LMI);
end

