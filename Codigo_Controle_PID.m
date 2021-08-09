%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VARIAVEIS ELETRICAS E MECANICAS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 0.373375;     % MASSA DO ROBO
r_d = 0.015;      % RAIO DA RODA DIREITA
r_e = 0.015;      % RAIO DA RODA ESQUERDA
d = 0.08;         % DIMENSÃO ROBO
K_d = 3.3*10^-2;  % CONSTANTE DE TORQUE MOTOR DIREITO
K_e = 3.3*10^-2;  % CONSTANTE DE TORQUE MOTOR ESQUERDO
J_d = 5.1*10^-6;  % MOMENTO DE INERCIA DO ACOPLAMENTO RODA/MOTOR DIREITO
J_e = 5.1*10^-6;  % MOMENTO DE INERCIA DO ACOPLAMENTO RODA/MOTOR ESQUERDO
B_d = 4.78*10^-6; % COEFICIENTE DE ATRITO DO MOTOR DIREITO
B_e = 4.78*10^-6; % COEFICIENTE DE ATRITO DO MOTOR ESQUERDO
p_d = 0.006;      % RESISTIVIDADE ELETRICA MOTOR DIREITO
p_e = 0.006;      % RESISTIVIDADE ELETRICA MOTOR ESQUERDO
J = 25.5*10^-6;   % MOMENTO DE INERCIA ROBO
B_t = 0.58*10^-6; % COEFICIENTE DE ATRITO COM O SOLO EM MOVIMENTO ROTACIONAL
B_l = 0.43;       % COEFICIENTE DE ATRITO COM O SOLO EM MOVIMENTO LINEAR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% MATRIZES %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J_m = [J_d   0;    % MATRIZ DE INERCIAS DOS ACOMPLAMENTOS RODAS/MOTORES
         0 J_e];
    
B_m = [B_d   0;    % MATRIZ DE COEFICIENTES DE ATRITO DOS MOTORES
         0 B_e];
    
J_r = [m  0;       % MATRIZ DE MASSA E INERCIA ROBO
       0  J];

B_r = [B_l   0;    % MATRIZ DE COEFICIENTES DE ATRITO SOLO/MOVIMENTOS
         0 B_t];

K_m = [K_d   0;    % MATRIZ DE CONSTANTES DE TORQUE DOS MOTORES
         0 K_e];    
    
p = [p_d    0;     % MATRIZ DE RESISTIVIDADE ELETRICA DOS MOTORES
       0  p_e];

w_T_u = [1/r_d  d/(2*r_d);  % MATRIZ DE CONVERSAO DE VELOCIDADES DO REFERENCIAL DO ROBO
         1/r_e -d/(2*r_e)]; % PARA O REFERENCIAL DOS ATUADORES
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% MATRIZES MODELO %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARTE ELETRICA
K_u = w_T_u' * p * K_m; 
% PARTE MECANICA
M_u = J_r + w_T_u' * J_m * w_T_u;
% RESTRICOES ATRITO/ACIONAMENTO
B_u = B_r + w_T_u' * [p * K_m * K_m + B_m] * w_T_u;   
% MATRIZ CONVERSAO SINAL DE CONTROLE ROT/LIN <-> TENSAO NOS MOTORES
T = -inv(-inv(M_u)*B_u)*inv(M_u)*K_u;                   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% MATRIZES MODELO CINEMATICO %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MODELO EM SS X_p = A*X + B*U
A = [-inv(M_u)*B_u zeros(2);
        eye(2)     zeros(2)];
    
B = [inv(M_u)*K_u ; zeros(2)] * inv(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% MODELOS %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = tf('s');

G_s = B(1,1) / (s^2 + B(1,1)*s) % MODELO LINEAR A PARTIR DE B
G_t = B(2,2) / (s^2 + B(2,2)*s) % MODELO ROTACIONAL A PARTIR DE B


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% CONTROLADOR LINEAR %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ta = 3.6;              % TEMPO DE ACOMODACAO DESEJADO
zeta = 1.1;            % ZETA PRA SISTEMA SUPERARMOTECIDO
Wn = 4 / (zeta*Ta);    % FREQ. NAT. DE OSC. DESEJADO

Kd_s = (2*zeta*Wn - B(1,1)) / B(1,1); % GANHO DERIVATIVO
Kp_s = (Wn*Wn)/B(1,1);                % GANHO PROPORCIONAL

C_s = Kp_s + Kd_s*s;                  % CONTROLADOR PD LINEAR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CONTROLADOR ANGULAR %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ta = 1.5;                       % TEMPO DE ACOMODACAO DESEJADO           
OV = 0.0001;                    % MAXIMO PICO DESEJADO
x = log(OV)^2;
zeta = sqrt ( x / (x + pi^2));  % ZETA PARA Ta e OV DESEJADOS
Wn = 4 / (zeta*Ta);             % FREQ. NAT. DE OSC. DESEJADO

Kd_t = (2*zeta*Wn - B(2,2)) / B(2,2); % GANHO DERIVATIVO
Kp_t = (Wn*Wn)/B(2,2);                % GANHO PROPORCIONAL

C_t = Kp_t + Kd_t*s;                  % CONTROLADOR PD ANGULAR





