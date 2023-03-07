clear;
close;

%% Definições
% Número de bits
num_b = 100000;
% Faixa de Eb/N0
Eb_N0_dB = 0:1:9;
% Faixa de Eb/N0 linearizada
Eb_N0_lin = 10 .^ (Eb_N0_dB/10);

%% Fonte
info = randi(2, 1, num_b) - 1;
% Adiciona um sufixo de zeros para finalizar as máquinas convolucionais no
% estado inicial
info(end-9:end) = zeros(1, 10);
    
%% Codificador

%%% Código Convolucional 1 com razão 1/2
convolutional_one_size = num_b * 2;
convolutional_one_info = zeros(1, convolutional_one_size);

for convolutional_one = 1
    convolutional_machine_one = containers.Map();

    % node = { edge_0 (next node, first output, second output), edge_1 (...) }
    convolutional_machine_one('00') = {{'00', 0, 0}, {'10', 1, 1}};
    convolutional_machine_one('10') = {{'01', 1, 0}, {'11', 0, 1}};
    convolutional_machine_one('11') = {{'01', 0, 1}, {'11', 1, 0}};
    convolutional_machine_one('01') = {{'00', 1, 1}, {'10', 0, 0}};

    current_state = '00';
    count = 1;
    for i = 1:num_b
        state_node = convolutional_machine_one(current_state);
        state_edge = state_node{info(i) + 1};
        current_state = state_edge{1};
        convolutional_one_info(count) = state_edge{2};
        count = count + 1;
        convolutional_one_info(count) = state_edge{3};
        count = count + 1;
    end
end

%%% Código Convolucional 2 com razão 1/3
convolutional_two_size = num_b * 3;
convolutional_two_info = zeros(1, convolutional_two_size);

for convolutional_two = 1
    convolutional_machine_two = containers.Map();

    % node = { edge_0 (next node, first output, second output), edge_1 (...) }
    convolutional_machine_two('00') = {{'00', 0, 0, 0}, {'10', 1, 1, 1}};
    convolutional_machine_two('10') = {{'01', 1, 1, 1}, {'11', 0, 0, 0}};
    convolutional_machine_two('11') = {{'01', 0, 0, 1}, {'11', 1, 1, 0}};
    convolutional_machine_two('01') = {{'00', 1, 1, 0}, {'10', 0, 0, 1}};

    current_state = '00';
    count = 1;
    for i = 1:num_b
        state_node = convolutional_machine_two(current_state);
        state_edge = state_node{info(i) + 1};
        current_state = state_edge{1};
        convolutional_two_info(count) = state_edge{2};
        count = count + 1;
        convolutional_two_info(count) = state_edge{3};
        count = count + 1;
        convolutional_two_info(count) = state_edge{4};
        count = count + 1;
    end
end

%%% Código Convolucional GSM Full Rate
convolutional_gsm_size = num_b * 2;
convolutional_gsm_info = zeros(1, convolutional_gsm_size);

for convolutional_gsm = 1
    convolutional_machine_gsm = containers.Map();

    % node = { edge_0 (next node, first output, second output), edge_1 (...) }
    convolutional_machine_gsm('a') = {{'a', 0, 0}, {'i', 1, 1}}; 
    convolutional_machine_gsm('b') = {{'a', 1, 1}, {'i', 0, 0}};
    convolutional_machine_gsm('c') = {{'b', 1, 1}, {'j', 0, 0}};
    convolutional_machine_gsm('d') = {{'b', 0, 0}, {'j', 1, 1}};
    convolutional_machine_gsm('e') = {{'c', 0, 0}, {'k', 1, 1}};
    convolutional_machine_gsm('f') = {{'c', 1, 1}, {'k', 0, 0}};
    convolutional_machine_gsm('g') = {{'d', 1, 1}, {'l', 0, 0}};
    convolutional_machine_gsm('h') = {{'d', 0, 0}, {'l', 1, 1}};
    convolutional_machine_gsm('i') = {{'e', 0, 1}, {'m', 1, 0}};
    convolutional_machine_gsm('j') = {{'e', 1, 0}, {'m', 0, 1}};
    convolutional_machine_gsm('k') = {{'f', 1, 0}, {'n', 0, 1}};
    convolutional_machine_gsm('l') = {{'f', 0, 1}, {'n', 1, 0}};
    convolutional_machine_gsm('m') = {{'g', 0, 1}, {'o', 1, 0}};
    convolutional_machine_gsm('n') = {{'g', 1, 0}, {'o', 0, 1}};
    convolutional_machine_gsm('o') = {{'h', 1, 0}, {'p', 0, 1}};
    convolutional_machine_gsm('p') = {{'h', 0, 1}, {'p', 1, 0}};

    current_state = 'a';
    count = 1;
    for i = 1:num_b
        state_node = convolutional_machine_gsm(current_state);
        state_edge = state_node{info(i) + 1};
        current_state = state_edge{1};
        convolutional_gsm_info(count) = state_edge{2};
        count = count + 1;
        convolutional_gsm_info(count) = state_edge{3};
        count = count + 1;
    end
end


%% Modulação

%%% BPSK
% Sem codificação
info_bpsk = complex(2*info-1, 0);
% Convolucional 1
convolutional_one_info_bpsk = complex(2*convolutional_one_info-1, 0);
% Convolucional 2
convolutional_two_info_bpsk = complex(2*convolutional_two_info-1, 0);
% Convolucional GSM
convolutional_gsm_info_bpsk = complex(2*convolutional_gsm_info-1, 0);

%%% 4-QAM
% Bits | Q          | I
% 00   | sqrt(2)/2  | sqrt(2)/2
% 01   | sqrt(2)/2  | -sqrt(2)/2
% 10   | -sqrt(2)/2 | sqrt(2)/2
% 11   | -sqrt(2)/2 | -sqrt(2)/2

% Sem codificação
qam_size = num_b / 2;
info_qam_I = zeros(1, qam_size);
info_qam_Q = zeros(1, qam_size);
for without_encoding = 1
    count = 1;
    for i = 1:2:(length(info)-1)
        if info(i) == 0
            if info(i+1) == 0 %% 00
                info_qam_I(count) = sqrt(2)/2;
                info_qam_Q(count) = sqrt(2)/2;
            else %% 01
                info_qam_I(count) = -1 * sqrt(2)/2;
                info_qam_Q(count) = sqrt(2)/2;
            end
        else 
            if info(i+1) == 0 %% 10
                info_qam_I(count) = sqrt(2)/2;
                info_qam_Q(count) = -1 * sqrt(2)/2;
            else %% 11
                info_qam_I(count) = -1 * sqrt(2)/2;
                info_qam_Q(count) = -1 * sqrt(2)/2;
            end
        end
        count = count + 1;
    end
end

% Convolucional 1
qam_size = convolutional_one_size / 2;
convolutional_one_info_qam_I = zeros(1, qam_size);
convolutional_one_info_qam_Q = zeros(1, qam_size);
for convolutional_one = 1
    count = 1;
    for i = 1:2:(convolutional_one_size-1)
        if convolutional_one_info(i) == 0
            if convolutional_one_info(i+1) == 0 %% 00
                convolutional_one_info_qam_I(count) = sqrt(2)/2;
                convolutional_one_info_qam_Q(count) = sqrt(2)/2;
            else %% 01
                convolutional_one_info_qam_I(count) = -1 * sqrt(2)/2;
                convolutional_one_info_qam_Q(count) = sqrt(2)/2;
            end
        else 
            if convolutional_one_info(i+1) == 0 %% 10
                convolutional_one_info_qam_I(count) = sqrt(2)/2;
                convolutional_one_info_qam_Q(count) = -1 * sqrt(2)/2;
            else %% 11
                convolutional_one_info_qam_I(count) = -1 * sqrt(2)/2;
                convolutional_one_info_qam_Q(count) = -1 * sqrt(2)/2;
            end
        end
        count = count + 1;
    end
end

% Convolucional 2
qam_size = convolutional_two_size / 2;
convolutional_two_info_qam_I = zeros(1, qam_size);
convolutional_two_info_qam_Q = zeros(1, qam_size);
for convolutional_one = 1
    count = 1;
    for i = 1:2:(convolutional_two_size-1)
        if convolutional_two_info(i) == 0
            if convolutional_two_info(i+1) == 0 %% 00
                convolutional_two_info_qam_I(count) = sqrt(2)/2;
                convolutional_two_info_qam_Q(count) = sqrt(2)/2;
            else %% 01
                convolutional_two_info_qam_I(count) = -1 * sqrt(2)/2;
                convolutional_two_info_qam_Q(count) = sqrt(2)/2;
            end
        else 
            if convolutional_two_info(i+1) == 0 %% 10
                convolutional_two_info_qam_I(count) = sqrt(2)/2;
                convolutional_two_info_qam_Q(count) = -1 * sqrt(2)/2;
            else %% 11
                convolutional_two_info_qam_I(count) = -1 * sqrt(2)/2;
                convolutional_two_info_qam_Q(count) = -1 * sqrt(2)/2;
            end
        end
        count = count + 1;
    end
end

% Convolucional GSM
qam_size = convolutional_gsm_size / 2;
convolutional_gsm_info_qam_I = zeros(1, qam_size);
convolutional_gsm_info_qam_Q = zeros(1, qam_size);
for convolutional_gsm = 1
    count = 1;
    for i = 1:2:(convolutional_gsm_size-1)
        if convolutional_gsm_info(i) == 0
            if convolutional_gsm_info(i+1) == 0 %% 00
                convolutional_gsm_info_qam_I(count) = sqrt(2)/2;
                convolutional_gsm_info_qam_Q(count) = sqrt(2)/2;
            else %% 01
                convolutional_gsm_info_qam_I(count) = -1 * sqrt(2)/2;
                convolutional_gsm_info_qam_Q(count) = sqrt(2)/2;
            end
        else 
            if convolutional_gsm_info(i+1) == 0 %% 10
                convolutional_gsm_info_qam_I(count) = sqrt(2)/2;
                convolutional_gsm_info_qam_Q(count) = -1 * sqrt(2)/2;
            else %% 11
                convolutional_gsm_info_qam_I(count) = -1 * sqrt(2)/2;
                convolutional_gsm_info_qam_Q(count) = -1 * sqrt(2)/2;
            end
        end
        count = count + 1;
    end
end


%% Receptor com BSPK

%%% Sem codificação
Eb = 1; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP); 
ber_bpsk_without_code = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, num_b), randn(1, num_b))*sqrt(0.5); 
   
    % Vetores recebido
    info_bpsk_with_noise = info_bpsk + n; 
        
    % Recupera a informação (sinal da parte real)
    real_info_with_noise = real(info_bpsk_with_noise);
    
    %%% Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    disp('BPSK sem codificação');
    disp(sum(info ~= demod));
                    
    % Contagem de erros e cálculo do BER
    ber_bpsk_without_code(i) = sum(info ~= demod) / num_b; 
end

%%% Convolucional 1
Eb = 2; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP); 
ber_bpsk_convolutional_one = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_one_size), randn(1, convolutional_one_size))*sqrt(0.5); 
   
    % Vetores recebido com ruído
    real_info_with_noise = real(convolutional_one_info_bpsk + n);
    
    %%% Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:2:length(demod)-1
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);

        % Armazena estados atuais
        viterbi_machine_00 = viterbi_machine('00');
        viterbi_machine_01 = viterbi_machine('01');
        viterbi_machine_10 = viterbi_machine('10');
        viterbi_machine_11 = viterbi_machine('11');
                
        % Análise de deslocamento para 00
        zero_zero_difference = Inf; % Valor máximo possível é 2
        zero_one_difference = Inf;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = viterbi_machine_00{2};
            if first_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = viterbi_machine_01{2};
            if first_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= Inf || zero_one_difference ~= Inf
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('00') = { [origin_state{1}, 0], difference };
        else
            viterbi_machine('00') = {};
        end
        
        % Análise de deslocamento para 01
        one_zero_difference = Inf;
        one_one_difference = Inf;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = viterbi_machine_10{2};
            if first_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = viterbi_machine_11{2};
            if first_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= Inf || one_one_difference ~= Inf
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('01') = { [origin_state{1}, 0], difference };
        else
            viterbi_machine('01') = {};
        end
        
        % Análise de deslocamento para 10
        zero_zero_difference = Inf;
        zero_one_difference = Inf;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = viterbi_machine_00{2};
            if first_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = viterbi_machine_01{2};
            if first_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= Inf || zero_one_difference ~= Inf
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('10') = { [origin_state{1}, 1], difference };
        else
            viterbi_machine('10') = {};
        end
        
        % Análise de deslocamento para 11        
        one_zero_difference = Inf;
        one_one_difference = Inf;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = viterbi_machine_10{2};
            if first_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = viterbi_machine_11{2};
            if first_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= Inf || one_one_difference ~= Inf
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('11') = { [origin_state{1}, 1], difference };
        else
            viterbi_machine('11') = {};
        end
    end
    
    % Seleciona o valor atual do estado 00
    zero_zero_state = viterbi_machine('00');

    decoded_sequence = zero_zero_state{1};
                
    disp('BPSK convolucional 1');
    disp(sum(info ~= decoded_sequence));
    
    % Contagem de erros e cálculo do BER
    ber_bpsk_convolutional_one(i) = sum(info ~= decoded_sequence) /  num_b; 
end

%%% Convolucional 2
Eb = 3; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP); 
ber_bpsk_convolutional_two = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_two_size), randn(1, convolutional_two_size))*sqrt(0.5); 
 
    % Recupera a informação com ruído (sinal da parte real)
    real_info_with_noise = real(convolutional_two_info_bpsk + n);
    
    %%% Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:3:length(demod)-2
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);
        third_bit = demod(x + 2);

        % Armazena estados atuais
        viterbi_machine_00 = viterbi_machine('00');
        viterbi_machine_01 = viterbi_machine('01');
        viterbi_machine_10 = viterbi_machine('10');
        viterbi_machine_11 = viterbi_machine('11');
                
        % Análise de deslocamento para 00
        zero_zero_difference = Inf; % Valor máximo possível é 3
        zero_one_difference = Inf;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = viterbi_machine_00{2};
            if first_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if third_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = viterbi_machine_01{2};
            if first_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
            if third_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= Inf || zero_one_difference ~= Inf
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('00') = { [origin_state{1}, 0], difference };
        else
            viterbi_machine('00') = {};
        end
        
        % Análise de deslocamento para 01
        one_zero_difference = Inf;
        one_one_difference = Inf;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = viterbi_machine_10{2};
            if first_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if third_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = viterbi_machine_11{2};
            if first_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if third_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= Inf || one_one_difference ~= Inf
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('01') = { [origin_state{1}, 0], difference };
        else
            viterbi_machine('01') = {};
        end
        
        % Análise de deslocamento para 10
        zero_zero_difference = Inf;
        zero_one_difference = Inf;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = viterbi_machine_00{2};
            if first_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if third_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = viterbi_machine_01{2};
            if first_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
            if third_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= Inf || zero_one_difference ~= Inf
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('10') = { [origin_state{1}, 1], difference };
        else
            viterbi_machine('10') = {};
        end
        
        % Análise de deslocamento para 11        
        one_zero_difference = Inf;
        one_one_difference = Inf;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = viterbi_machine_10{2};
            if first_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if third_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = viterbi_machine_11{2};
            if first_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if third_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= Inf || one_one_difference ~= Inf
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('11') = { [origin_state{1}, 1], difference };
        else
            viterbi_machine('11') = {};
        end
    end
    
    % Seleciona o valor atual do estado 00
    zero_zero_state = viterbi_machine('00');
    decoded_sequence = zero_zero_state{1};
            
    disp('BPSK convolucional 2');
    disp(sum(info ~= decoded_sequence));
    
    % Contagem de erros e cálculo do BER
    ber_bpsk_convolutional_two(i) = sum(info ~= decoded_sequence) / num_b; 
end

%%% Convolucional GSM
Eb = 2; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_bpsk_convolutional_gsm = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_gsm_size), randn(1, convolutional_gsm_size))*sqrt(0.5); 
 
    % Recupera a informação com ruído (sinal da parte real)
    real_info_with_noise = real(convolutional_gsm_info_bpsk + n); 
    
    %%% Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('a') = {[], 0};
    viterbi_machine('b') = {};
    viterbi_machine('c') = {};
    viterbi_machine('d') = {};
    viterbi_machine('e') = {};
    viterbi_machine('f') = {};
    viterbi_machine('g') = {};
    viterbi_machine('h') = {};
    viterbi_machine('i') = {};
    viterbi_machine('j') = {};
    viterbi_machine('k') = {};
    viterbi_machine('l') = {};
    viterbi_machine('m') = {};
    viterbi_machine('n') = {};
    viterbi_machine('o') = {};
    viterbi_machine('p') = {};
    
    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:2:length(demod)-1
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);

        % Armazena estados atuais
        previous_viterbi_machine = containers.Map();
        for copy_viterbi_machine = 1
            previous_viterbi_machine('a') = viterbi_machine('a');
            previous_viterbi_machine('b') = viterbi_machine('b');
            previous_viterbi_machine('c') = viterbi_machine('c');
            previous_viterbi_machine('d') = viterbi_machine('d');
            previous_viterbi_machine('e') = viterbi_machine('e');
            previous_viterbi_machine('f') = viterbi_machine('f');
            previous_viterbi_machine('g') = viterbi_machine('g');
            previous_viterbi_machine('h') = viterbi_machine('h');
            previous_viterbi_machine('i') = viterbi_machine('i');
            previous_viterbi_machine('j') = viterbi_machine('j');
            previous_viterbi_machine('k') = viterbi_machine('k');
            previous_viterbi_machine('l') = viterbi_machine('l');
            previous_viterbi_machine('m') = viterbi_machine('m');
            previous_viterbi_machine('n') = viterbi_machine('n');
            previous_viterbi_machine('o') = viterbi_machine('o');
            previous_viterbi_machine('p') = viterbi_machine('p');
        end
        
        % Análise de deslocamento para a
        for state_a = 1
            a_difference = Inf;
            b_difference = Inf;
            a_state = previous_viterbi_machine('a');
            b_state = previous_viterbi_machine('b');
            if ~isempty(a_state)
                a_difference = a_state{2};
                if first_bit ~= 0
                    a_difference = a_difference + 1;
                end
                if second_bit ~= 0
                    a_difference = a_difference + 1;
                end
            end
            if ~isempty(b_state)
                b_difference = b_state{2};
                if first_bit ~= 1
                    b_difference = b_difference + 1;
                end
                if second_bit ~= 1
                    b_difference = b_difference + 1;
                end
            end
            if a_difference ~= Inf || b_difference ~= Inf
                if a_difference <= b_difference
                    new_bits = [a_state{1}, 0];
                    new_difference = a_difference;
                else
                    new_bits = [b_state{1}, 0];
                    new_difference = b_difference;
                end
                viterbi_machine('a') = { new_bits, new_difference };
            else
                viterbi_machine('a') = {};
            end
        end
        
        % Análise de deslocamento para b
        for state_b = 1
            c_difference = Inf;
            d_difference = Inf;
            c_state = previous_viterbi_machine('c');
            d_state = previous_viterbi_machine('d');
            if ~isempty(c_state)
                c_difference = c_state{2};
                if first_bit ~= 1
                    c_difference = c_difference + 1;
                end
                if second_bit ~= 1
                    c_difference = c_difference + 1;
                end
            end
            if ~isempty(d_state)
                d_difference = d_state{2};
                if first_bit ~= 0
                    d_difference = d_difference + 1;
                end
                if second_bit ~= 0
                    d_difference = d_difference + 1;
                end
            end
            if c_difference ~= Inf || d_difference ~= Inf
                if c_difference <= d_difference
                    new_bits = [c_state{1}, 0];
                    new_difference = c_difference;
                else
                    new_bits = [d_state{1}, 0];
                    new_difference = d_difference;
                end
                viterbi_machine('b') = { new_bits, new_difference };
            else
                viterbi_machine('b') = {};
            end
        end
        
        % Análise de deslocamento para c
        for state_c = 1
            e_difference = Inf;
            f_difference = Inf;
            e_state = previous_viterbi_machine('e');
            f_state = previous_viterbi_machine('f');
            if ~isempty(e_state)
                e_difference = e_state{2};
                if first_bit ~= 0
                    e_difference = e_difference + 1;
                end
                if second_bit ~= 0
                    e_difference = e_difference + 1;
                end
            end
            if ~isempty(f_state)
                f_difference = f_state{2};
                if first_bit ~= 1
                    f_difference = f_difference + 1;
                end
                if second_bit ~= 1
                    f_difference = f_difference + 1;
                end
            end
            if e_difference ~= Inf || f_difference ~= Inf
                if e_difference <= f_difference
                    new_bits = [e_state{1}, 0];
                    new_difference = e_difference;
                else
                    new_bits = [f_state{1}, 0];
                    new_difference = f_difference;
                end
                viterbi_machine('c') = { new_bits, new_difference };
            else
                viterbi_machine('c') = {};
            end
        end
        
        % Análise de deslocamento para d
        for state_d = 1
            g_difference = Inf;
            h_difference = Inf;
            g_state = previous_viterbi_machine('g');
            h_state = previous_viterbi_machine('h');
            if ~isempty(g_state)
                g_difference = g_state{2};
                if first_bit ~= 1
                    g_difference = g_difference + 1;
                end
                if second_bit ~= 1
                    g_difference = g_difference + 1;
                end
            end
            if ~isempty(h_state)
                h_difference = h_state{2};
                if first_bit ~= 0
                    h_difference = h_difference + 1;
                end
                if second_bit ~= 0
                    h_difference = h_difference + 1;
                end
            end
            if g_difference ~= Inf || h_difference ~= Inf
                if g_difference <= h_difference
                    new_bits = [g_state{1}, 0];
                    new_difference = g_difference;
                else
                    new_bits = [h_state{1}, 0];
                    new_difference = h_difference;
                end
                viterbi_machine('d') = { new_bits, new_difference };
            else
                viterbi_machine('d') = {};
            end
        end
        
        % Análise de deslocamento para e
        for state_e = 1
            i_difference = Inf;
            j_difference = Inf;
            i_state = previous_viterbi_machine('i');
            j_state = previous_viterbi_machine('j');
            if ~isempty(i_state)
                i_difference = i_state{2};
                if first_bit ~= 0
                    i_difference = i_difference + 1;
                end
                if second_bit ~= 1
                    i_difference = i_difference + 1;
                end
            end
            if ~isempty(j_state)
                j_difference = j_state{2};
                if first_bit ~= 1
                    j_difference = j_difference + 1;
                end
                if second_bit ~= 0
                    j_difference = j_difference + 1;
                end
            end
            if i_difference ~= Inf || j_difference ~= Inf
                if i_difference <= j_difference
                    new_bits = [i_state{1}, 0];
                    new_difference = i_difference;
                else
                    new_bits = [j_state{1}, 0];
                    new_difference = j_difference;
                end
                viterbi_machine('e') = { new_bits, new_difference };
            else
                viterbi_machine('e') = {};
            end
        end
        
        % Análise de deslocamento para f
        for state_f = 1
            k_difference = Inf;
            l_difference = Inf;
            k_state = previous_viterbi_machine('k');
            l_state = previous_viterbi_machine('l');
            if ~isempty(k_state)
                k_difference = k_state{2};
                if first_bit ~= 1
                    k_difference = k_difference + 1;
                end
                if second_bit ~= 0
                    k_difference = k_difference + 1;
                end
            end
            if ~isempty(l_state)
                l_difference = l_state{2};
                if first_bit ~= 0
                    l_difference = l_difference + 1;
                end
                if second_bit ~= 1
                    l_difference = l_difference + 1;
                end
            end
            if k_difference ~= Inf || l_difference ~= Inf
                if k_difference <= l_difference
                    new_bits = [k_state{1}, 0];
                    new_difference = k_difference;
                else
                    new_bits = [l_state{1}, 0];
                    new_difference = l_difference;
                end
                viterbi_machine('f') = { new_bits, new_difference };
            else
                viterbi_machine('f') = {};
            end
        end
        
        % Análise de deslocamento para g
        for state_g = 1
            m_difference = Inf;
            n_difference = Inf;
            m_state = previous_viterbi_machine('m');
            n_state = previous_viterbi_machine('n');
            if ~isempty(m_state)
                m_difference = m_state{2};
                if first_bit ~= 0
                    m_difference = m_difference + 1;
                end
                if second_bit ~= 1
                    m_difference = m_difference + 1;
                end
            end
            if ~isempty(n_state)
                n_difference = n_state{2};
                if first_bit ~= 1
                    n_difference = n_difference + 1;
                end
                if second_bit ~= 0
                    n_difference = n_difference + 1;
                end
            end
            if m_difference ~= Inf || n_difference ~= Inf
                if m_difference <= n_difference
                    new_bits = [m_state{1}, 0];
                    new_difference = m_difference;
                else
                    new_bits = [n_state{1}, 0];
                    new_difference = n_difference;
                end
                viterbi_machine('g') = { new_bits, new_difference };
            else
                viterbi_machine('g') = {};
            end
        end
        
        % Análise de deslocamento para h
        for state_h = 1
            o_difference = Inf;
            p_difference = Inf;
            o_state = previous_viterbi_machine('o');
            p_state = previous_viterbi_machine('p');
            if ~isempty(o_state)
                o_difference = o_state{2};
                if first_bit ~= 1
                    o_difference = o_difference + 1;
                end
                if second_bit ~= 0
                    o_difference = o_difference + 1;
                end
            end
            if ~isempty(p_state)
                p_difference = p_state{2};
                if first_bit ~= 0
                    p_difference = p_difference + 1;
                end
                if second_bit ~= 1
                    p_difference = p_difference + 1;
                end
            end
            if o_difference ~= Inf || p_difference ~= Inf
                if o_difference <= p_difference
                    new_bits = [o_state{1}, 0];
                    new_difference = o_difference;
                else
                    new_bits = [p_state{1}, 0];
                    new_difference = p_difference;
                end
                viterbi_machine('h') = { new_bits, new_difference };
            else
                viterbi_machine('h') = {};
            end
        end
        
        % Análise de deslocamento para i
        for state_i = 1
            a_difference = Inf;
            b_difference = Inf;
            a_state = previous_viterbi_machine('a');
            b_state = previous_viterbi_machine('b');
            if ~isempty(a_state)
                a_difference = a_state{2};
                if first_bit ~= 1
                    a_difference = a_difference + 1;
                end
                if second_bit ~= 1
                    a_difference = a_difference + 1;
                end
            end
            if ~isempty(b_state)
                b_difference = b_state{2};
                if first_bit ~= 0
                    b_difference = b_difference + 1;
                end
                if second_bit ~= 0
                    b_difference = b_difference + 1;
                end
            end
            if a_difference ~= Inf || b_difference ~= Inf
                if a_difference <= b_difference
                    new_bits = [a_state{1}, 1];
                    new_difference = a_difference;
                else
                    new_bits = [b_state{1}, 1];
                    new_difference = b_difference;
                end
                viterbi_machine('i') = { new_bits, new_difference };
            else
                viterbi_machine('i') = {};
            end
        end
        
        % Análise de deslocamento para j
        for state_j = 1
            c_difference = Inf;
            d_difference = Inf;
            c_state = previous_viterbi_machine('c');
            d_state = previous_viterbi_machine('d');
            if ~isempty(c_state)
                c_difference = c_state{2};
                if first_bit ~= 0
                    c_difference = c_difference + 1;
                end
                if second_bit ~= 0
                    c_difference = c_difference + 1;
                end
            end
            if ~isempty(d_state)
                d_difference = d_state{2};
                if first_bit ~= 1
                    d_difference = d_difference + 1;
                end
                if second_bit ~= 1
                    d_difference = d_difference + 1;
                end
            end
            if c_difference ~= Inf || d_difference ~= Inf
                if c_difference <= d_difference
                    new_bits = [c_state{1}, 1];
                    new_difference = c_difference;
                else
                    new_bits = [d_state{1}, 1];
                    new_difference = d_difference;
                end
                viterbi_machine('j') = { new_bits, new_difference };
            else
                viterbi_machine('j') = {};
            end
        end
        
        % Análise de deslocamento para k
        for state_k = 1
            e_difference = Inf;
            f_difference = Inf;
            e_state = previous_viterbi_machine('e');
            f_state = previous_viterbi_machine('f');
            if ~isempty(e_state)
                e_difference = e_state{2};
                if first_bit ~= 1
                    e_difference = e_difference + 1;
                end
                if second_bit ~= 1
                    e_difference = e_difference + 1;
                end
            end
            if ~isempty(f_state)
                f_difference = f_state{2};
                if first_bit ~= 0
                    f_difference = f_difference + 1;
                end
                if second_bit ~= 0
                    f_difference = f_difference + 1;
                end
            end
            if e_difference ~= Inf || f_difference ~= Inf
                if e_difference <= f_difference
                    new_bits = [e_state{1}, 1];
                    new_difference = e_difference;
                else
                    new_bits = [f_state{1}, 1];
                    new_difference = f_difference;
                end
                viterbi_machine('k') = { new_bits, new_difference };
            else
                viterbi_machine('k') = {};
            end
        end
        
        % Análise de deslocamento para l
        for state_l = 1
            g_difference = Inf;
            h_difference = Inf;
            g_state = previous_viterbi_machine('g');
            h_state = previous_viterbi_machine('h');
            if ~isempty(g_state)
                g_difference = g_state{2};
                if first_bit ~= 0
                    g_difference = g_difference + 1;
                end
                if second_bit ~= 0
                    g_difference = g_difference + 1;
                end
            end
            if ~isempty(h_state)
                h_difference = h_state{2};
                if first_bit ~= 1
                    h_difference = h_difference + 1;
                end
                if second_bit ~= 1
                    h_difference = h_difference + 1;
                end
            end
            if g_difference ~= Inf || h_difference ~= Inf
                if g_difference <= h_difference
                    new_bits = [g_state{1}, 1];
                    new_difference = g_difference;
                else
                    new_bits = [h_state{1}, 1];
                    new_difference = h_difference;
                end
                viterbi_machine('l') = { new_bits, new_difference };
            else
                viterbi_machine('l') = {};
            end
        end
        
        % Análise de deslocamento para m
        for state_m = 1
            i_difference = Inf;
            j_difference = Inf;
            i_state = previous_viterbi_machine('i');
            j_state = previous_viterbi_machine('j');
            if ~isempty(i_state)
                i_difference = i_state{2};
                if first_bit ~= 1
                    i_difference = i_difference + 1;
                end
                if second_bit ~= 0
                    i_difference = i_difference + 1;
                end
            end
            if ~isempty(j_state)
                j_difference = j_state{2};
                if first_bit ~= 0
                    j_difference = j_difference + 1;
                end
                if second_bit ~= 1
                    j_difference = j_difference + 1;
                end
            end
            if i_difference ~= Inf || j_difference ~= Inf
                if i_difference <= j_difference
                    new_bits = [i_state{1}, 1];
                    new_difference = i_difference;
                else
                    new_bits = [j_state{1}, 1];
                    new_difference = j_difference;
                end
                viterbi_machine('m') = { new_bits, new_difference };
            else
                viterbi_machine('m') = {};
            end
        end
        
        % Análise de deslocamento para n
        for state_n = 1
            k_difference = Inf;
            l_difference = Inf;
            k_state = previous_viterbi_machine('k');
            l_state = previous_viterbi_machine('l');
            if ~isempty(k_state)
                k_difference = k_state{2};
                if first_bit ~= 0
                    k_difference = k_difference + 1;
                end
                if second_bit ~= 1
                    k_difference = k_difference + 1;
                end
            end
            if ~isempty(l_state)
                l_difference = l_state{2};
                if first_bit ~= 1
                    l_difference = l_difference + 1;
                end
                if second_bit ~= 0
                    l_difference = l_difference + 1;
                end
            end
            if k_difference ~= Inf || l_difference ~= Inf
                if k_difference <= l_difference
                    new_bits = [k_state{1}, 1];
                    new_difference = k_difference;
                else
                    new_bits = [l_state{1}, 1];
                    new_difference = l_difference;
                end
                viterbi_machine('n') = { new_bits, new_difference };
            else
                viterbi_machine('n') = {};
            end
        end
        
        % Análise de deslocamento para o
        for state_o = 1
            m_difference = Inf;
            n_difference = Inf;
            m_state = previous_viterbi_machine('m');
            n_state = previous_viterbi_machine('n');
            if ~isempty(m_state)
                m_difference = m_state{2};
                if first_bit ~= 1
                    m_difference = m_difference + 1;
                end
                if second_bit ~= 0
                    m_difference = m_difference + 1;
                end
            end
            if ~isempty(n_state)
                n_difference = n_state{2};
                if first_bit ~= 0
                    n_difference = n_difference + 1;
                end
                if second_bit ~= 1
                    n_difference = n_difference + 1;
                end
            end
            if m_difference ~= Inf || n_difference ~= Inf
                if m_difference <= n_difference
                    new_bits = [m_state{1}, 1];
                    new_difference = m_difference;
                else
                    new_bits = [n_state{1}, 1];
                    new_difference = n_difference;
                end
                viterbi_machine('o') = { new_bits, new_difference };
            else
                viterbi_machine('o') = {};
            end
        end
        
        % Análise de deslocamento para p
        for state_p = 1
            o_difference = Inf;
            p_difference = Inf;
            o_state = previous_viterbi_machine('o');
            p_state = previous_viterbi_machine('p');
            if ~isempty(o_state)
                o_difference = o_state{2};
                if first_bit ~= 0
                    o_difference = o_difference + 1;
                end
                if second_bit ~= 1
                    o_difference = o_difference + 1;
                end
            end
            if ~isempty(p_state)
                p_difference = p_state{2};
                if first_bit ~= 1
                    p_difference = p_difference + 1;
                end
                if second_bit ~= 0
                    p_difference = p_difference + 1;
                end
            end
            if o_difference ~= Inf || p_difference ~= Inf
                if o_difference <= p_difference
                    new_bits = [o_state{1}, 1];
                    new_difference = o_difference;
                else
                    new_bits = [p_state{1}, 1];
                    new_difference = p_difference;
                end
                viterbi_machine('p') = { new_bits, new_difference };
            else
                viterbi_machine('p') = {};
            end
        end
    end
    
    % Seleciona o valor atual do estado 00
    initial_state = viterbi_machine('a');
    decoded_sequence = initial_state{1};
    
    disp('BPSK GSM');
    disp(sum(info ~= decoded_sequence));
    
    % Contagem de erros e cálculo do BER
    ber_bpsk_convolutional_gsm(i) = sum(info ~= decoded_sequence) / num_b; 
end

%% Receptor com 4QAM
 
%%% Sem codificação
Eb = 1.5; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_qam_without_code = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, num_b / 2), randn(1, num_b / 2))*sqrt(0.5); 
   
    % Vetores recebido
    I_with_noise = real(info_qam_I + n);
    Q_with_noise = real(info_qam_Q + n);
     
    %%% Demodulação            
    demod = zeros(1, length(I_with_noise) * 2);
    count = 1;
    for x = 1:length(I_with_noise)
        if I_with_noise(x) >= 0
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 0;
            else
                demod(count) = 1;
                demod(count + 1) = 0;
            end
        else
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 1;
            else
                demod(count) = 1;
                demod(count + 1) = 1;
            end
        end
        count = count + 2;
    end
    
    disp('qam sem códificação');
    disp(sum(info ~= demod));   
    
    % Contagem de erros e cálculo do BER
    ber_qam_without_code(i) = sum(info ~= demod) / num_b; 
end

%%% Convolucional 1
Eb = 3; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_qam_convolutional_one = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_one_size / 2), randn(1, convolutional_one_size / 2))*sqrt(0.5); 
   
    % Vetores recebido
    I_with_noise = real(convolutional_one_info_qam_I + n); 
    Q_with_noise = real(convolutional_one_info_qam_Q + n); 
     
    %%% Demodulação            
    demod = zeros(1, length(I_with_noise) * 2);
    count = 1;
    for x = 1:length(I_with_noise)
        if I_with_noise(x) >= 0
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 0;
            else
                demod(count) = 1;
                demod(count + 1) = 0;
            end
        else
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 1;
            else
                demod(count) = 1;
                demod(count + 1) = 1;
            end
        end
        count = count + 2;
    end
    
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:2:length(demod)-1
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);

        % Armazena estados atuais
        viterbi_machine_00 = viterbi_machine('00');
        viterbi_machine_01 = viterbi_machine('01');
        viterbi_machine_10 = viterbi_machine('10');
        viterbi_machine_11 = viterbi_machine('11');
                
        % Análise de deslocamento para 00
        zero_zero_difference = Inf; % Valor máximo possível é 2
        zero_one_difference = Inf;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = viterbi_machine_00{2};
            if first_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = viterbi_machine_01{2};
            if first_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= Inf || zero_one_difference ~= Inf
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('00') = { [origin_state{1}, 0], difference };
        else
            viterbi_machine('00') = {};
        end
        
        % Análise de deslocamento para 01
        one_zero_difference = Inf;
        one_one_difference = Inf;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = viterbi_machine_10{2};
            if first_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = viterbi_machine_11{2};
            if first_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= Inf || one_one_difference ~= Inf
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('01') = { [origin_state{1}, 0], difference };
        else
            viterbi_machine('01') = {};
        end
        
        % Análise de deslocamento para 10
        zero_zero_difference = Inf;
        zero_one_difference = Inf;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = viterbi_machine_00{2};
            if first_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = viterbi_machine_01{2};
            if first_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= Inf || zero_one_difference ~= Inf
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('10') = { [origin_state{1}, 1], difference };
        else
            viterbi_machine('10') = {};
        end
        
        % Análise de deslocamento para 11        
        one_zero_difference = Inf;
        one_one_difference = Inf;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = viterbi_machine_10{2};
            if first_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = viterbi_machine_11{2};
            if first_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= Inf || one_one_difference ~= Inf
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('11') = { [origin_state{1}, 1], difference };
        else
            viterbi_machine('11') = {};
        end
    end
    
    % Seleciona o valor atual do estado 00
    zero_zero_state = viterbi_machine('00');
    decoded_sequence = zero_zero_state{1};
       
    disp('qam convolucional 1');
    disp(sum(info ~= decoded_sequence));   
    
    % Contagem de erros e cálculo do BER
    ber_qam_convolutional_one(i) = sum(info ~= decoded_sequence) / num_b; 
end

%%% Convolucional 2
Eb = 4.5; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_qam_convolutional_two = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_two_size / 2), randn(1, convolutional_two_size / 2))*sqrt(0.5); 
   
    % Vetores recebido
    I_with_noise = real(convolutional_two_info_qam_I + n); 
    Q_with_noise = real(convolutional_two_info_qam_Q + n); 
     
    %%% Demodulação            
    demod = zeros(1, length(I_with_noise) * 2);
    count = 1;
    for x = 1:length(I_with_noise)
        if I_with_noise(x) >= 0
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 0;
            else
                demod(count) = 1;
                demod(count + 1) = 0;
            end
        else
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 1;
            else
                demod(count) = 1;
                demod(count + 1) = 1;
            end
        end
        count = count + 2;
    end
        
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:3:length(demod)-2
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);
        third_bit = demod(x + 2);

        % Armazena estados atuais
        viterbi_machine_00 = viterbi_machine('00');
        viterbi_machine_01 = viterbi_machine('01');
        viterbi_machine_10 = viterbi_machine('10');
        viterbi_machine_11 = viterbi_machine('11');
                
        % Análise de deslocamento para 00
        zero_zero_difference = Inf; % Valor máximo possível é 3
        zero_one_difference = Inf;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = viterbi_machine_00{2};
            if first_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if third_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = viterbi_machine_01{2};
            if first_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
            if third_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= Inf || zero_one_difference ~= Inf
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('00') = { [origin_state{1}, 0], difference };
        else
            viterbi_machine('00') = {};
        end
        
        % Análise de deslocamento para 01
        one_zero_difference = Inf;
        one_one_difference = Inf;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = viterbi_machine_10{2};
            if first_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if third_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = viterbi_machine_11{2};
            if first_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if third_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= Inf || one_one_difference ~= Inf
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('01') = { [origin_state{1}, 0], difference };
        else
            viterbi_machine('01') = {};
        end
        
        % Análise de deslocamento para 10
        zero_zero_difference = Inf;
        zero_one_difference = Inf;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = viterbi_machine_00{2};
            if first_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if third_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = viterbi_machine_01{2};
            if first_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
            if third_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= Inf || zero_one_difference ~= Inf
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('10') = { [origin_state{1}, 1], difference };
        else
            viterbi_machine('10') = {};
        end
        
        % Análise de deslocamento para 11        
        one_zero_difference = Inf;
        one_one_difference = Inf;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = viterbi_machine_10{2};
            if first_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if third_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = viterbi_machine_11{2};
            if first_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if third_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= Inf || one_one_difference ~= Inf
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('11') = { [origin_state{1}, 1], difference };
        else
            viterbi_machine('11') = {};
        end
    end
    
    % Seleciona o valor atual do estado 00
    zero_zero_state = viterbi_machine('00');
    decoded_sequence = zero_zero_state{1};
       
    disp('qam convolucional 2');
    disp(sum(info ~= decoded_sequence));  
    
    % Contagem de erros e cálculo do BER
    ber_qam_convolutional_two(i) = sum(info ~= decoded_sequence) / num_b; 
end

%%% Convolucional GSM
Eb = 3; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_qam_convolutional_gsm = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_gsm_size / 2), randn(1, convolutional_gsm_size / 2))*sqrt(0.5); 
   
    % Vetores recebido
    I_with_noise = real(convolutional_gsm_info_qam_I + n); 
    Q_with_noise = real(convolutional_gsm_info_qam_Q + n); 
     
    %%% Demodulação            
    demod = zeros(1, length(I_with_noise) * 2);
    count = 1;
    for x = 1:length(I_with_noise)
        if I_with_noise(x) >= 0
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 0;
            else
                demod(count) = 1;
                demod(count + 1) = 0;
            end
        else
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 1;
            else
                demod(count) = 1;
                demod(count + 1) = 1;
            end
        end
        count = count + 2;
    end
        
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('a') = {[], 0};
    viterbi_machine('b') = {};
    viterbi_machine('c') = {};
    viterbi_machine('d') = {};
    viterbi_machine('e') = {};
    viterbi_machine('f') = {};
    viterbi_machine('g') = {};
    viterbi_machine('h') = {};
    viterbi_machine('i') = {};
    viterbi_machine('j') = {};
    viterbi_machine('k') = {};
    viterbi_machine('l') = {};
    viterbi_machine('m') = {};
    viterbi_machine('n') = {};
    viterbi_machine('o') = {};
    viterbi_machine('p') = {};
    
    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:2:length(demod)-1
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);

        % Armazena estados atuais
        previous_viterbi_machine = containers.Map();
        for copy_viterbi_machine = 1
            previous_viterbi_machine('a') = viterbi_machine('a');
            previous_viterbi_machine('b') = viterbi_machine('b');
            previous_viterbi_machine('c') = viterbi_machine('c');
            previous_viterbi_machine('d') = viterbi_machine('d');
            previous_viterbi_machine('e') = viterbi_machine('e');
            previous_viterbi_machine('f') = viterbi_machine('f');
            previous_viterbi_machine('g') = viterbi_machine('g');
            previous_viterbi_machine('h') = viterbi_machine('h');
            previous_viterbi_machine('i') = viterbi_machine('i');
            previous_viterbi_machine('j') = viterbi_machine('j');
            previous_viterbi_machine('k') = viterbi_machine('k');
            previous_viterbi_machine('l') = viterbi_machine('l');
            previous_viterbi_machine('m') = viterbi_machine('m');
            previous_viterbi_machine('n') = viterbi_machine('n');
            previous_viterbi_machine('o') = viterbi_machine('o');
            previous_viterbi_machine('p') = viterbi_machine('p');
        end
        
        % Análise de deslocamento para a
        for state_a = 1
            a_difference = Inf;
            b_difference = Inf;
            a_state = previous_viterbi_machine('a');
            b_state = previous_viterbi_machine('b');
            if ~isempty(a_state)
                a_difference = a_state{2};
                if first_bit ~= 0
                    a_difference = a_difference + 1;
                end
                if second_bit ~= 0
                    a_difference = a_difference + 1;
                end
            end
            if ~isempty(b_state)
                b_difference = b_state{2};
                if first_bit ~= 1
                    b_difference = b_difference + 1;
                end
                if second_bit ~= 1
                    b_difference = b_difference + 1;
                end
            end
            if a_difference ~= Inf || b_difference ~= Inf
                if a_difference <= b_difference
                    new_bits = [a_state{1}, 0];
                    new_difference = a_difference;
                else
                    new_bits = [b_state{1}, 0];
                    new_difference = b_difference;
                end
                viterbi_machine('a') = { new_bits, new_difference };
            else
                viterbi_machine('a') = {};
            end
        end
        
        % Análise de deslocamento para b
        for state_b = 1
            c_difference = Inf;
            d_difference = Inf;
            c_state = previous_viterbi_machine('c');
            d_state = previous_viterbi_machine('d');
            if ~isempty(c_state)
                c_difference = c_state{2};
                if first_bit ~= 1
                    c_difference = c_difference + 1;
                end
                if second_bit ~= 1
                    c_difference = c_difference + 1;
                end
            end
            if ~isempty(d_state)
                d_difference = d_state{2};
                if first_bit ~= 0
                    d_difference = d_difference + 1;
                end
                if second_bit ~= 0
                    d_difference = d_difference + 1;
                end
            end
            if c_difference ~= Inf || d_difference ~= Inf
                if c_difference <= d_difference
                    new_bits = [c_state{1}, 0];
                    new_difference = c_difference;
                else
                    new_bits = [d_state{1}, 0];
                    new_difference = d_difference;
                end
                viterbi_machine('b') = { new_bits, new_difference };
            else
                viterbi_machine('b') = {};
            end
        end
        
        % Análise de deslocamento para c
        for state_c = 1
            e_difference = Inf;
            f_difference = Inf;
            e_state = previous_viterbi_machine('e');
            f_state = previous_viterbi_machine('f');
            if ~isempty(e_state)
                e_difference = e_state{2};
                if first_bit ~= 0
                    e_difference = e_difference + 1;
                end
                if second_bit ~= 0
                    e_difference = e_difference + 1;
                end
            end
            if ~isempty(f_state)
                f_difference = f_state{2};
                if first_bit ~= 1
                    f_difference = f_difference + 1;
                end
                if second_bit ~= 1
                    f_difference = f_difference + 1;
                end
            end
            if e_difference ~= Inf || f_difference ~= Inf
                if e_difference <= f_difference
                    new_bits = [e_state{1}, 0];
                    new_difference = e_difference;
                else
                    new_bits = [f_state{1}, 0];
                    new_difference = f_difference;
                end
                viterbi_machine('c') = { new_bits, new_difference };
            else
                viterbi_machine('c') = {};
            end
        end
        
        % Análise de deslocamento para d
        for state_d = 1
            g_difference = Inf;
            h_difference = Inf;
            g_state = previous_viterbi_machine('g');
            h_state = previous_viterbi_machine('h');
            if ~isempty(g_state)
                g_difference = g_state{2};
                if first_bit ~= 1
                    g_difference = g_difference + 1;
                end
                if second_bit ~= 1
                    g_difference = g_difference + 1;
                end
            end
            if ~isempty(h_state)
                h_difference = h_state{2};
                if first_bit ~= 0
                    h_difference = h_difference + 1;
                end
                if second_bit ~= 0
                    h_difference = h_difference + 1;
                end
            end
            if g_difference ~= Inf || h_difference ~= Inf
                if g_difference <= h_difference
                    new_bits = [g_state{1}, 0];
                    new_difference = g_difference;
                else
                    new_bits = [h_state{1}, 0];
                    new_difference = h_difference;
                end
                viterbi_machine('d') = { new_bits, new_difference };
            else
                viterbi_machine('d') = {};
            end
        end
        
        % Análise de deslocamento para e
        for state_e = 1
            i_difference = Inf;
            j_difference = Inf;
            i_state = previous_viterbi_machine('i');
            j_state = previous_viterbi_machine('j');
            if ~isempty(i_state)
                i_difference = i_state{2};
                if first_bit ~= 0
                    i_difference = i_difference + 1;
                end
                if second_bit ~= 1
                    i_difference = i_difference + 1;
                end
            end
            if ~isempty(j_state)
                j_difference = j_state{2};
                if first_bit ~= 1
                    j_difference = j_difference + 1;
                end
                if second_bit ~= 0
                    j_difference = j_difference + 1;
                end
            end
            if i_difference ~= Inf || j_difference ~= Inf
                if i_difference <= j_difference
                    new_bits = [i_state{1}, 0];
                    new_difference = i_difference;
                else
                    new_bits = [j_state{1}, 0];
                    new_difference = j_difference;
                end
                viterbi_machine('e') = { new_bits, new_difference };
            else
                viterbi_machine('e') = {};
            end
        end
        
        % Análise de deslocamento para f
        for state_f = 1
            k_difference = Inf;
            l_difference = Inf;
            k_state = previous_viterbi_machine('k');
            l_state = previous_viterbi_machine('l');
            if ~isempty(k_state)
                k_difference = k_state{2};
                if first_bit ~= 1
                    k_difference = k_difference + 1;
                end
                if second_bit ~= 0
                    k_difference = k_difference + 1;
                end
            end
            if ~isempty(l_state)
                l_difference = l_state{2};
                if first_bit ~= 0
                    l_difference = l_difference + 1;
                end
                if second_bit ~= 1
                    l_difference = l_difference + 1;
                end
            end
            if k_difference ~= Inf || l_difference ~= Inf
                if k_difference <= l_difference
                    new_bits = [k_state{1}, 0];
                    new_difference = k_difference;
                else
                    new_bits = [l_state{1}, 0];
                    new_difference = l_difference;
                end
                viterbi_machine('f') = { new_bits, new_difference };
            else
                viterbi_machine('f') = {};
            end
        end
        
        % Análise de deslocamento para g
        for state_g = 1
            m_difference = Inf;
            n_difference = Inf;
            m_state = previous_viterbi_machine('m');
            n_state = previous_viterbi_machine('n');
            if ~isempty(m_state)
                m_difference = m_state{2};
                if first_bit ~= 0
                    m_difference = m_difference + 1;
                end
                if second_bit ~= 1
                    m_difference = m_difference + 1;
                end
            end
            if ~isempty(n_state)
                n_difference = n_state{2};
                if first_bit ~= 1
                    n_difference = n_difference + 1;
                end
                if second_bit ~= 0
                    n_difference = n_difference + 1;
                end
            end
            if m_difference ~= Inf || n_difference ~= Inf
                if m_difference <= n_difference
                    new_bits = [m_state{1}, 0];
                    new_difference = m_difference;
                else
                    new_bits = [n_state{1}, 0];
                    new_difference = n_difference;
                end
                viterbi_machine('g') = { new_bits, new_difference };
            else
                viterbi_machine('g') = {};
            end
        end
        
        % Análise de deslocamento para h
        for state_h = 1
            o_difference = Inf;
            p_difference = Inf;
            o_state = previous_viterbi_machine('o');
            p_state = previous_viterbi_machine('p');
            if ~isempty(o_state)
                o_difference = o_state{2};
                if first_bit ~= 1
                    o_difference = o_difference + 1;
                end
                if second_bit ~= 0
                    o_difference = o_difference + 1;
                end
            end
            if ~isempty(p_state)
                p_difference = p_state{2};
                if first_bit ~= 0
                    p_difference = p_difference + 1;
                end
                if second_bit ~= 1
                    p_difference = p_difference + 1;
                end
            end
            if o_difference ~= Inf || p_difference ~= Inf
                if o_difference <= p_difference
                    new_bits = [o_state{1}, 0];
                    new_difference = o_difference;
                else
                    new_bits = [p_state{1}, 0];
                    new_difference = p_difference;
                end
                viterbi_machine('h') = { new_bits, new_difference };
            else
                viterbi_machine('h') = {};
            end
        end
        
        % Análise de deslocamento para i
        for state_i = 1
            a_difference = Inf;
            b_difference = Inf;
            a_state = previous_viterbi_machine('a');
            b_state = previous_viterbi_machine('b');
            if ~isempty(a_state)
                a_difference = a_state{2};
                if first_bit ~= 1
                    a_difference = a_difference + 1;
                end
                if second_bit ~= 1
                    a_difference = a_difference + 1;
                end
            end
            if ~isempty(b_state)
                b_difference = b_state{2};
                if first_bit ~= 0
                    b_difference = b_difference + 1;
                end
                if second_bit ~= 0
                    b_difference = b_difference + 1;
                end
            end
            if a_difference ~= Inf || b_difference ~= Inf
                if a_difference <= b_difference
                    new_bits = [a_state{1}, 1];
                    new_difference = a_difference;
                else
                    new_bits = [b_state{1}, 1];
                    new_difference = b_difference;
                end
                viterbi_machine('i') = { new_bits, new_difference };
            else
                viterbi_machine('i') = {};
            end
        end
        
        % Análise de deslocamento para j
        for state_j = 1
            c_difference = Inf;
            d_difference = Inf;
            c_state = previous_viterbi_machine('c');
            d_state = previous_viterbi_machine('d');
            if ~isempty(c_state)
                c_difference = c_state{2};
                if first_bit ~= 0
                    c_difference = c_difference + 1;
                end
                if second_bit ~= 0
                    c_difference = c_difference + 1;
                end
            end
            if ~isempty(d_state)
                d_difference = d_state{2};
                if first_bit ~= 1
                    d_difference = d_difference + 1;
                end
                if second_bit ~= 1
                    d_difference = d_difference + 1;
                end
            end
            if c_difference ~= Inf || d_difference ~= Inf
                if c_difference <= d_difference
                    new_bits = [c_state{1}, 1];
                    new_difference = c_difference;
                else
                    new_bits = [d_state{1}, 1];
                    new_difference = d_difference;
                end
                viterbi_machine('j') = { new_bits, new_difference };
            else
                viterbi_machine('j') = {};
            end
        end
        
        % Análise de deslocamento para k
        for state_k = 1
            e_difference = Inf;
            f_difference = Inf;
            e_state = previous_viterbi_machine('e');
            f_state = previous_viterbi_machine('f');
            if ~isempty(e_state)
                e_difference = e_state{2};
                if first_bit ~= 1
                    e_difference = e_difference + 1;
                end
                if second_bit ~= 1
                    e_difference = e_difference + 1;
                end
            end
            if ~isempty(f_state)
                f_difference = f_state{2};
                if first_bit ~= 0
                    f_difference = f_difference + 1;
                end
                if second_bit ~= 0
                    f_difference = f_difference + 1;
                end
            end
            if e_difference ~= Inf || f_difference ~= Inf
                if e_difference <= f_difference
                    new_bits = [e_state{1}, 1];
                    new_difference = e_difference;
                else
                    new_bits = [f_state{1}, 1];
                    new_difference = f_difference;
                end
                viterbi_machine('k') = { new_bits, new_difference };
            else
                viterbi_machine('k') = {};
            end
        end
        
        % Análise de deslocamento para l
        for state_l = 1
            g_difference = Inf;
            h_difference = Inf;
            g_state = previous_viterbi_machine('g');
            h_state = previous_viterbi_machine('h');
            if ~isempty(g_state)
                g_difference = g_state{2};
                if first_bit ~= 0
                    g_difference = g_difference + 1;
                end
                if second_bit ~= 0
                    g_difference = g_difference + 1;
                end
            end
            if ~isempty(h_state)
                h_difference = h_state{2};
                if first_bit ~= 1
                    h_difference = h_difference + 1;
                end
                if second_bit ~= 1
                    h_difference = h_difference + 1;
                end
            end
            if g_difference ~= Inf || h_difference ~= Inf
                if g_difference <= h_difference
                    new_bits = [g_state{1}, 1];
                    new_difference = g_difference;
                else
                    new_bits = [h_state{1}, 1];
                    new_difference = h_difference;
                end
                viterbi_machine('l') = { new_bits, new_difference };
            else
                viterbi_machine('l') = {};
            end
        end
        
        % Análise de deslocamento para m
        for state_m = 1
            i_difference = Inf;
            j_difference = Inf;
            i_state = previous_viterbi_machine('i');
            j_state = previous_viterbi_machine('j');
            if ~isempty(i_state)
                i_difference = i_state{2};
                if first_bit ~= 1
                    i_difference = i_difference + 1;
                end
                if second_bit ~= 0
                    i_difference = i_difference + 1;
                end
            end
            if ~isempty(j_state)
                j_difference = j_state{2};
                if first_bit ~= 0
                    j_difference = j_difference + 1;
                end
                if second_bit ~= 1
                    j_difference = j_difference + 1;
                end
            end
            if i_difference ~= Inf || j_difference ~= Inf
                if i_difference <= j_difference
                    new_bits = [i_state{1}, 1];
                    new_difference = i_difference;
                else
                    new_bits = [j_state{1}, 1];
                    new_difference = j_difference;
                end
                viterbi_machine('m') = { new_bits, new_difference };
            else
                viterbi_machine('m') = {};
            end
        end
        
        % Análise de deslocamento para n
        for state_n = 1
            k_difference = Inf;
            l_difference = Inf;
            k_state = previous_viterbi_machine('k');
            l_state = previous_viterbi_machine('l');
            if ~isempty(k_state)
                k_difference = k_state{2};
                if first_bit ~= 0
                    k_difference = k_difference + 1;
                end
                if second_bit ~= 1
                    k_difference = k_difference + 1;
                end
            end
            if ~isempty(l_state)
                l_difference = l_state{2};
                if first_bit ~= 1
                    l_difference = l_difference + 1;
                end
                if second_bit ~= 0
                    l_difference = l_difference + 1;
                end
            end
            if k_difference ~= Inf || l_difference ~= Inf
                if k_difference <= l_difference
                    new_bits = [k_state{1}, 1];
                    new_difference = k_difference;
                else
                    new_bits = [l_state{1}, 1];
                    new_difference = l_difference;
                end
                viterbi_machine('n') = { new_bits, new_difference };
            else
                viterbi_machine('n') = {};
            end
        end
        
        % Análise de deslocamento para o
        for state_o = 1
            m_difference = Inf;
            n_difference = Inf;
            m_state = previous_viterbi_machine('m');
            n_state = previous_viterbi_machine('n');
            if ~isempty(m_state)
                m_difference = m_state{2};
                if first_bit ~= 1
                    m_difference = m_difference + 1;
                end
                if second_bit ~= 0
                    m_difference = m_difference + 1;
                end
            end
            if ~isempty(n_state)
                n_difference = n_state{2};
                if first_bit ~= 0
                    n_difference = n_difference + 1;
                end
                if second_bit ~= 1
                    n_difference = n_difference + 1;
                end
            end
            if m_difference ~= Inf || n_difference ~= Inf
                if m_difference <= n_difference
                    new_bits = [m_state{1}, 1];
                    new_difference = m_difference;
                else
                    new_bits = [n_state{1}, 1];
                    new_difference = n_difference;
                end
                viterbi_machine('o') = { new_bits, new_difference };
            else
                viterbi_machine('o') = {};
            end
        end
        
        % Análise de deslocamento para p
        for state_p = 1
            o_difference = Inf;
            p_difference = Inf;
            o_state = previous_viterbi_machine('o');
            p_state = previous_viterbi_machine('p');
            if ~isempty(o_state)
                o_difference = o_state{2};
                if first_bit ~= 0
                    o_difference = o_difference + 1;
                end
                if second_bit ~= 1
                    o_difference = o_difference + 1;
                end
            end
            if ~isempty(p_state)
                p_difference = p_state{2};
                if first_bit ~= 1
                    p_difference = p_difference + 1;
                end
                if second_bit ~= 0
                    p_difference = p_difference + 1;
                end
            end
            if o_difference ~= Inf || p_difference ~= Inf
                if o_difference <= p_difference
                    new_bits = [o_state{1}, 1];
                    new_difference = o_difference;
                else
                    new_bits = [p_state{1}, 1];
                    new_difference = p_difference;
                end
                viterbi_machine('p') = { new_bits, new_difference };
            else
                viterbi_machine('p') = {};
            end
        end
    end
    
    % Seleciona o valor atual do estado 00
    initial_state = viterbi_machine('a');
    decoded_sequence = initial_state{1};
    
    disp('qam GSM');
    disp(sum(info ~= decoded_sequence));
    
    % Contagem de erros e cálculo do BER
    ber_qam_convolutional_gsm(i) = sum(info ~= decoded_sequence) / num_b; 
end

%% Avaliação resultado

% BER teórico para comparação
semilogy(Eb_N0_dB, ber_qam_convolutional_gsm, 'x', Eb_N0_dB, ber_bpsk_convolutional_gsm, 'x', Eb_N0_dB, ber_qam_without_code, 'x', Eb_N0_dB, ber_qam_convolutional_one, 'x', Eb_N0_dB, ber_qam_convolutional_two, 'x', Eb_N0_dB, ber_bpsk_without_code, 'x', Eb_N0_dB, ber_bpsk_convolutional_one, 'x', Eb_N0_dB, ber_bpsk_convolutional_two, 'x', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('4QAM GSM', 'BPSK GSM', '4QAM sem codificação', '4QAM com código convolucional de razão 1/2', '4QAM com código convolucional de razão 1/3', 'BPSK sem codificação', 'BPSK com código convolucional de razão 1/2', 'BPSK com código convolucional de razão 1/3');
