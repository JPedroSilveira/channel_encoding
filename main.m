clear;
close;

%% Definições
% Número de bits
num_b = 100000;
% faixa de Eb/N0
Eb_N0_dB = 0:1:9;
% faixa de Eb/N0 linearizada
Eb_N0_lin = 10 .^ (Eb_N0_dB/10);

%% Fonte
info = randi(2, 1, num_b) - 1;

%% Codificador

%%% Viterbi
convolutional_machine = containers.Map();

% node = { edge_0 (next node, first output, second output), edge_1 }
convolutional_machine('00') = {{'00', 0, 0}, {'10', 1, 1}};
convolutional_machine('10') = {{'01', 1, 0}, {'11', 0, 1}};
convolutional_machine('11') = {{'01', 0, 1}, {'11', 1, 0}};
convolutional_machine('01') = {{'00', 1, 1}, {'10', 0, 0}};

convolutional_end = containers.Map();
convolutional_end('00') = [0, 0, 0, 0, 0, 0];
convolutional_end('10') = [0, 1, 0, 1, 1, 1];
convolutional_end('11') = [0, 1, 1, 1, 0, 0];
convolutional_end('01') = [1, 1, 0, 0, 0, 0];

current_state = '00';
viterbi_info_size = num_b * 2 + 6;
viterbi_info = zeros(1, viterbi_info_size);
count = 1;

for i = 1:num_b
    state_node = convolutional_machine(current_state);
    state_edge = state_node{info(i) + 1};
    current_state = state_edge{1};
    viterbi_info(count) = state_edge{2};
    count = count + 1;
    viterbi_info(count) = state_edge{3};
    count = count + 1;
end

viterbi_suffix = convolutional_end(current_state);
for i = 1:length(viterbi_suffix)
    viterbi_info(count) = viterbi_suffix(i);
    count = count + 1;
end

%% Modulação

%%% BPSK
viterbi_info_bpsk = complex(2*viterbi_info-1, 0);

%%% 4-QAM
% Bits | Q          | I
% 00   | sqrt(2)/2  | sqrt(2)/2
% 01   | sqrt(2)/2  | -sqrt(2)/2
% 10   | -sqrt(2)/2 | sqrt(2)/2
% 11   | -sqrt(2)/2 | -sqrt(2)/2

qam_size = viterbi_info_size / 2;
I = zeros(1, qam_size);
Q = zeros(1, qam_size);
count = 1;
for b = 1:2:(viterbi_info_size-1)
    if viterbi_info(b) == 0
        if viterbi_info(b+1) == 0 %% 00
            I(count) = sqrt(2)/2;
            Q(count) = sqrt(2)/2;
        else %% 01
            I(count) = -1 * sqrt(2)/2;
            Q(count) = sqrt(2)/2;
        end
    else 
        if viterbi_info(b+1) == 0 %% 10
            I(count) = sqrt(2)/2;
            Q(count) = -1 * sqrt(2)/2;
        else %% 11
            I(count) = -1 * sqrt(2)/2;
            Q(count) = -1 * sqrt(2)/2;
        end
    end
    count = count + 1;
end
viterbi_info_4qam = I + Q;

%% Receptor com BSPK

% Energia por bit para a modulação BPSK utilizada
Eb = 1; 
% Vetor de potências do ruído
NP = Eb ./ (Eb_N0_lin); 
% Vetor de amplitudes do ruído
NA = sqrt(NP); 

%%% Viterbi

for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, viterbi_info_size), randn(1, viterbi_info_size))*sqrt(0.5); 
   
    % Vetores recebido
    viterbi_info_bpsk_with_noise = viterbi_info_bpsk + n; 
        
    % Recupera a informação (sinal da parte real)
    real_info_with_noise = real(viterbi_info_bpsk_with_noise);
    
    % Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    % Decodificação
    current_state = '00';
    viterbi_state = containers.Map();
    % node = { edge_0 (next node, first output, second output), edge_1 }
    viterbi_state('00') = {{'00', 0, 0}, {'10', 1, 1}};
    viterbi_state('10') = {{'01', 1, 0}, {'11', 0, 1}};
    viterbi_state('11') = {{'01', 0, 1}, {'11', 1, 0}};
    viterbi_state('01') = {{'00', 1, 1}, {'10', 0, 0}};
    
    disp(demod);
                
    % Contagem de erros e cálculo do BER
    ber(i) = sum(bits ~= demod) / viterbi_info_size; 
end

%% Avaliação resultado

%% Prof
% bits aleatórios modulados em BPSK (parte real em 1 e -1)
bits = complex(2*randi(2, 1, num_b)-3, 0); 

% faixa de Eb/N0
Eb_N0_dB = 0:1:9; 

% faixa de Eb/N0 linearizada
Eb_N0_lin = 10 .^ (Eb_N0_dB/10);

% pré-alocação do vetor de BER
ber = zeros(size(Eb_N0_lin)); 

% energia por bit para a modulação BPSK utilizada
Eb = 1; 

% vetor de potências do ruído
NP = Eb ./ (Eb_N0_lin); 

% vetor de amplitudes do ruído
NA = sqrt(NP); 
    
for i = 1:length(Eb_N0_lin)
    % vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, num_b), randn(1, num_b))*sqrt(0.5); 
   
    % vetor recebido
    r = bits + n; 
        
    % recupera a informação (sinal da parte real)
    demod = sign(real(r)); 
            
    % contagem de erros e cálculo do BER
    ber(i) = sum(bits ~= demod) / num_b; 
end

%BER teórico para comparação
ber_theoretical = 0.5*erfc(sqrt(2*Eb_N0_lin)/sqrt(2)); 

semilogy(Eb_N0_dB, ber, 'x', Eb_N0_dB, ber_theoretical, 'r', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('Simulado','Teórico');
