clear;
close;

%% Definições
% Número de bits
num_b = 1000000;
% Faixa de Eb/N0
Eb_N0_dB = 0:1:12;
% Faixa de Eb/N0 linearizada
Eb_N0_lin = 10 .^ (Eb_N0_dB/10);
% Viterbi profundidade
tb = 16;

%% Fonte
info = randi(2, 1, num_b) - 1;
% Adiciona um sufixo de zeros para finalizar as máquinas convolucionais no
% estado inicial
info(end-9:end) = zeros(1, 10);
    
%% Codificador

%%% Código Convolucional 1 com razão 1/2
convolutional_one_size = num_b * 2;
trellis_convolutional_one = poly2trellis(3, [7 5]);
convolutional_one_info = convenc(info, trellis_convolutional_one);

%%% Código Convolucional 2 com razão 1/3
convolutional_two_size = num_b * 3;
trellis_convolutional_two = poly2trellis(3, [4 6 7]);
convolutional_two_info = convenc(info, trellis_convolutional_two);

%%% Código Convolucional GSM Full Rate
convolutional_gsm_size = num_b * 2;
trellis_convolutional_gsm = poly2trellis(5, [23 27]);
convolutional_gsm_info = convenc(info, trellis_convolutional_gsm);

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
    decoded_sequence = vitdec(demod,trellis_convolutional_one,tb,'trunc','hard');
                
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
    decoded_sequence = vitdec(demod,trellis_convolutional_two,tb,'trunc','hard');
            
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
    decoded_sequence = vitdec(demod,trellis_convolutional_gsm,tb,'trunc','hard');
    
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
    decoded_sequence = vitdec(demod,trellis_convolutional_one,tb,'trunc','hard');
       
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
    decoded_sequence = vitdec(demod,trellis_convolutional_two,tb,'trunc','hard');
       
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
    decoded_sequence = vitdec(demod,trellis_convolutional_gsm,tb,'trunc','hard');
        
    disp('qam GSM');
    disp(sum(info ~= decoded_sequence));
    
    % Contagem de erros e cálculo do BER
    ber_qam_convolutional_gsm(i) = sum(info ~= decoded_sequence) / num_b; 
end

%% Avaliação resultado

% BER teórico para comparação
semilogy(Eb_N0_dB, ber_qam_convolutional_gsm, 'r', Eb_N0_dB, ber_bpsk_convolutional_gsm, 'g', Eb_N0_dB, ber_qam_without_code, 'b', Eb_N0_dB, ber_qam_convolutional_one, 'c', Eb_N0_dB, ber_qam_convolutional_two, 'm', Eb_N0_dB, ber_bpsk_without_code, 'y', Eb_N0_dB, ber_bpsk_convolutional_one, 'b', Eb_N0_dB, ber_bpsk_convolutional_two, 'r', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('4QAM GSM', 'BPSK GSM', '4QAM sem codificação', '4QAM com código convolucional de razão 1/2', '4QAM com código convolucional de razão 1/3', 'BPSK sem codificação', 'BPSK com código convolucional de razão 1/2', 'BPSK com código convolucional de razão 1/3');
