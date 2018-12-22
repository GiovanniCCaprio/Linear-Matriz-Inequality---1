%% Lista B - Ex 1
clc,clear all,close all

% Dados

A{1} = [0 1 0; 0 0 1; -1 -3 -2];
A{2} = [0 1 0; 0 0 1; -1 -22 -2];
H = [1 -1;-1  1]; % lá em baixo é multiplicado por gamma...

% Configurações

order = size(A{1},1);
output.cpusec_m = clock;
val = [];
gamma = 0.001;
k = 1;
p = 0;

while p >= 0 
    % Criação das Variaveis
    P1 = sdpvar(order,order,'symmetric'); % Lyapunov > 0
    P2 = sdpvar(order,order,'symmetric'); % Lyapunov > 0
    G1 = sdpvar(order,order); % Não será usado? G e F
    G2 = sdpvar(order,order);
    F = sdpvar(order,order);
 
    % Variação a cada novo gamma dos h
    h1 = gamma*H(:,1);
    h2 = gamma*H(:,2);

    % Criação das LMIs
    LMIs = []; % novas LMIs a cada iteração
    LMIs = [LMIs, P1 >= 0];
    LMIs = [LMIs, P2 >= 0];

    L1 = A{1}'*P1 + P1*A{1} + h1(1,1)*P1 + h1(2,1)*P2;
    LMIs = [LMIs, L1 <= 0];
    L12 = A{1}'*P1 + P1*A{1} + h2(1,1)*P1 + h2(2,1)*P2;
    LMIs = [LMIs, L12 <= 0];

    L2 = A{2}'*P2 + P2*A{2} + h1(1,1)*P1 + h1(2,1)*P2;
    LMIs = [LMIs, L2 <= 0];
    L22 = A{2}'*P2 + P2*A{2} + h2(1,1)*P1 + h2(2,1)*P2;
    LMIs = [LMIs, L22 <= 0];

    L3 = A{1}'*P2 + P2*A{1} + A{2}'*P1 + P1*A{2} + 2*(h1(1,1)*P1 + h1(2,1)*P2);
    LMIs = [LMIs, L3 <= 0];
    L32 = A{1}'*P2 + P2*A{1} + A{2}'*P1 + P1*A{2} + 2*(h2(1,1)*P1 + h2(2,1)*P2);
    LMIs = [LMIs, L3 <= 0];

    % Número de Variáveis
    output.V = size(getvariables(LMIs),2);

    % Solução a LMIs
    sol = solvesdp(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));

    % Tempo para solução
    output.cpusec_s = sol.solvertime;

    % Resíduo
    p = min(checkset(LMIs));
    val{k} = p; % salva os residuos
    k = k + 1; % nova iteração
    gamma = 0.001 + k*0.01; % gamma novo
end

display('O gamma máximo é:')
gamma_max = gamma - 0.01  % pega o gamma anterior