clear variables; 
close all; 
clc
%% definition of some constant parameters 
% the channels
channels = {'C3', 'C4', 'CZ', 'F3', 'F4', 'F7', 'F8', 'FP1', 'FP2', 'FZ', 'O1', ...
    'O2', 'P3', 'P4', 'PZ', 'T3', 'T4', 'T5', 'T6'};
% sampling frequency
fs = 500;
% in order to iterate through all the subjects
subject_indexes = [1, 2, 3, 5, 7, 8];
base_subject = 'Subject0';
subject_count = 0;
max_N = 0;
% additional structure to compute the signals' mean at zones
zones_keys = {'central', 'frontal', 'simmetrical anterior frontal', 'occipital', 'parietal', 'temporal' };
zones_channels = {{'CZ', 'C3', 'C4'}, {'F3', 'F4', 'F7', 'F8','FZ'}, {'FP1', 'FP2'}, {'O1', 'O2'}, {'P3', 'P4', 'PZ'},{'T3', 'T4', 'T5', 'T6'}};
zones_map = containers.Map(zones_keys, zones_channels);
% bands in Hz
bands = struct('alpha', [8 12], 'beta', [13 30], 'theta', [4 7]);
% to store the CSA values for each subject for each channel use a matrix of
% cells
resting_PSD_list = cell(length(subject_indexes), 1);
math_PSD_list = cell(length(subject_indexes), 1);
freq_resting = [];
freq_math = [];
max_N_resting = 0;
max_N_math = 0;
% to have the frequency at resting and computing math
shift_ratio = 0.5;  % 50% overlap
nfft        = 2048;
noverlap    = fs * shift_ratio;

%% computation of the EEG graph and PSD for each channel for all subjects

for i = subject_indexes
    subject_index = int2str(i);
    current_subject = append(base_subject, subject_index);
    subject_count = subject_count + 1;
    
    % to distinguish the two phases under analysis
    % 1 = resting
    % 2 = math
    for j = 1 : 2
        current_file = append(current_subject, '_', int2str(j), '.mat' );
        % loading the current file
        load(current_file)
        disp(current_file);
        % container for the EEGs registered from the various channels
        EEG_list = [];
        % container for the PSD list
        tmp_PSD_list = cell(length(channels),1);
        max_N = 0;
        % iterate on each channels to build the data
        for k = 1 : length(channels)
            channel_name = channels{k};
            EEG = eval(channel_name);
            N = length(EEG);
            if N > max_N
                % define the time axis
                t = (1 : N)/fs;
                max_N = N;
            elseif N < max_N
                zeros_to_add = max_N - N;
                EEG = [ EEG; zeros(zeros_to_add, 1)];
            end
            % for the current EEG calculate the CSA 
            % [current_PSD, current_f] = pwelch(EEG, hamming(fs), noverlap, nfft, fs);
            % current_PSD = current_PSD(current_f < 30);
            % tmp_PSD_list is an array with as many rows as the number of
            % channels
            % for each channel k it contains a matrix whose rows are the
            % number of temporal segments (time) on which the PSD was
            % calculated, while the columns are the number of frequencies
            % (under 30 Hz) over which the PSD was calculated
            [~, freq] = pwelch(EEG, hamming(fs), noverlap, nfft, fs);
            freq = freq(freq < 30);
            tmp_PSD_list{k} = computeCSA(EEG, j);
            if j == 1
                freq_resting = freq;
                if max_N > max_N_resting
                    max_N_resting = max_N;
                end
            else
                freq_math = freq;
                if max_N > max_N_math
                    max_N_math = max_N;
                end
            end
            % adding the channel measurement to the EEG list for the subject
            EEG_list = [EEG_list EEG];

        end
        % computing some parameters depending on the phase
        if j == 1
            current_state = 'resting';
            resting_PSD_list{subject_count} = tmp_PSD_list;
        else 
            current_state = 'performing mathematical task';
            math_PSD_list{subject_count} = tmp_PSD_list;
        end
        % commented to make the code faster at the end remove comments
        figure
        % stacked plot of the EEGs
        channelLabels = "Channel " + channels;
        plotTitle = append('Subject ', subject_index, ' : ', current_state);
        stackedplot (t, EEG_list, 'b', 'DisplayLabels', channelLabels, 'Title', plotTitle);
        xlabel('Time [s]');
    end
end
disp('End computation part 1');

%% parameters for next phase
subject_count = 0;
time_resting = computeTime(1, max_N_resting);
time_math = computeTime(2, max_N_math);
pie_colors = [1, 0.31, 0.53;  
              1, 0.85, 0.46;  
              0.47, 0.89, 0.78]; 

%% computation of total PSD and alpha, beta and theta powers per zone
for i = subject_indexes
    subject_count = subject_count + 1;
    
    % for every zone, define a figure and create tile layout and set title
    % for every phase then 
    % - set total_PSD, alpha, beta, theta powers to zero
    % - consider all the channels in that zone
    %       add contribution to total_PSD, alpha, beta, theta
    % - show PSD and compute data for piechart in two consecutive tiles
    for zone = zones_keys
        figure
        % create a tiled layout of 2x2 plots
        tile = tiledlayout(2,2);
        overall_caption = append('Subject ', int2str(i), ' ', zone, ' zone');
        disp(overall_caption);
        sgtitle(overall_caption, 'FontSize', 14, 'FontWeight', 'bold');
        channels_in_zone = zones_map(zone{1});
        for j = 1 : 2
            if j == 1
                phase = 'resting';
                time = time_resting;
                frequency = freq_resting;
            else
                phase = 'math';
                time = time_math;
                frequency = freq_math;
            end
            % calculating the portion of the matrix for the various bands
            alpha_idx = frequency >= bands.alpha(1) & frequency <= bands.alpha(2);
            beta_idx = frequency >= bands.beta(1) & frequency <= bands.beta(2);
            theta_idx = frequency >= bands.theta(1) & frequency <= bands.theta(2);
            total_PSD = zeros(size(time, 1), size(freq_math, 1));
            total_alpha = 0;
            total_beta = 0;
            total_theta = 0;
            for ch = 1 : length(channels_in_zone)
                ch_index = find(strcmp(channels, channel_name));     
                if j == 1
                    CSA = resting_PSD_list{subject_count}{ch};
                else
                    CSA = math_PSD_list{subject_count}{ch};
                end
                total_alpha = total_alpha + mean(sum(CSA(:, alpha_idx), 2));
                total_beta = total_beta + mean(sum(CSA(:, beta_idx), 2));
                total_theta = total_theta + mean(sum(CSA(:, theta_idx), 2));
                total_PSD = total_PSD + CSA;
            end
            average_PSD = total_PSD / length(channels_in_zone);
            tot = total_alpha + total_beta + total_theta;
            alpha_percentage = total_alpha / tot * 100;
            beta_percentage = total_beta / tot * 100;
            theta_percentage = total_theta / tot * 100;
            
            nexttile
            % waterfall
            [wtime, wfreq] = meshgrid(time, frequency);
            waterfall(wtime', wfreq', total_PSD);
            colormap('turbo');  % more colorful
            shading interp;     % shadow
            string = sprintf('PSD during %s', phase);
            title(string);
            xlabel('Time [m]');
            ylabel('Frequency [Hz]');
            zlabel('\muV^{2}/Hz');
            grid on;

            nexttile

            %pie chart
            labels = {sprintf('Alpha: %.1f%%', alpha_percentage), ...
                      sprintf('Beta: %.1f%%', beta_percentage), ...
                      sprintf('Theta: %.1f%%', theta_percentage)};
            h = pie([alpha_percentage, beta_percentage, theta_percentage], labels);
            
            % display total power below the chart and add title
            tot_string = sprintf('Total Power: %.2f \\muV^{2}/Hz', tot);
            text(-1, -1.5, tot_string, 'HorizontalAlignment', 'center', 'FontSize', 12);
            string = append('Waves power during ', phase);
            title(string,'FontSize', 12);

            % assign colors to slices
            for k = 1:2:length(h)  % for every other element in the pie chart
                h(k).FaceColor = pie_colors((k+1)/2, :); % apply the corresponding color
            end

        end
    end
end

%% computing MSC
load('chanlocs.mat')
band = {'alpha', 'beta', 'theta'};
MSC_tot_math =zeros(length(channels), length(channels),length(band)); % creating a matrix of zeros to hold the sum oh the MSC of all the subjects in the three bands to then compute the mean value during arithmetic tasks
MSC_tot_resting = zeros(length(channels), length(channels),length(band)); % creating a matrix of zeros to hold the sum oh the MSC of all the subjects in the three bands to then compute the mean value during resting state
for i = subject_indexes
    subject_index = int2str(i);
    current_subject = append(base_subject, subject_index);
    
    
    % to distinguish the two phases under analysis
    % 1 = resting
    % 2 = math
    EEG_list_math = [];
    EEG_list_resting=[];
    for j = 1:2
    
        current_file = append(current_subject, '_', int2str(j), '.mat' );
        % loading the current file
        load(current_file)        
        max_N = 0;
        for k = 1 : length(channels)
            channel_name = channels{k};
            EEG = eval(channel_name);
            N = length(EEG);
            if N > max_N
                % define the time axis
                t = (1 : N)/fs;
                max_N = N;
            elseif N < max_N
                zeros_to_add = max_N - N;
                EEG = [ EEG; zeros(zeros_to_add, 1)];
            end
            
            if j==1
               EEG_list_resting = [EEG_list_resting EEG];
            else 
                EEG_list_math = [EEG_list_math EEG];
            end

        end
    end


    % computing the MSC for each subject for each band in both conditions
    % (resting state and arithmetic tasks) and the rapresenting it with imagesc
    % and topoplot. Then computing the average value of MSC of all subject for
    % each band in both condition and the rapresnting it.
    fs = 500;
    L_resting = 3 * 500; % window to compute the coherence for the EEG during resting
    L_math = 1 * 500; % window compute the coherence for the EEG during arithmetics tasks
    nfft = 2048;

    bands = struct('alpha', [8 12], 'beta', [13 30], 'theta', [4 7]);
    MSC_results_resting = zeros(size(EEG_list_resting,2), size(EEG_list_resting,2)); %creating a matrix of zerosto hold the results of the MSC during resting state for each subject
    band_names = fieldnames(bands);
    num_bands = length(band_names);

    for ind = 1:num_bands % calculating the MSC for each band and creating a MSC results matrix for each band
        figure;
        band_name = band_names{ind};
        frequency_range = bands.(band_name);
        freq_min = frequency_range(1);
        freq_max = frequency_range(2);

        % calculating the MSC for each signals during resting 
        for k = 1:size(EEG_list_resting, 2)
            for l = k:size(EEG_list_resting, 2)
                [MSC_resting, f] = mscohere(EEG_list_resting(:, k), EEG_list_resting(:, l), hamming(L_resting), round(0.5 * L_resting), nfft, fs);
                freq_indexes = f >= freq_min & f <= freq_max;
                filtered_MSC_resting = MSC_resting(freq_indexes); % extracting the MSC values of the specific band
                mean_MSC_resting = mean(filtered_MSC_resting); %computing the mean value for the band
                MSC_results_resting(k,l) = mean_MSC_resting; %storing the mean value in the results matrix (19x19)
                MSC_results_resting(l,k) = MSC_results_resting(k,l); %the results matrix is symmetrical
            end

        end
        MSC_tot_resting(:, :, ind) = MSC_tot_resting(:, :, ind) + MSC_results_resting; %summing the value of MSC during resting state of this band of all subjects to then compute the mean value

        MSC_results_math = zeros(size(EEG_list_math,2), size(EEG_list_math,2)); %creating a matrix of zeros to hold the results of the MSC during arithmetics tasks
        
        % calculating the MSC for each signals during arithmetics tasks
        for m = 1:size(EEG_list_math, 2)
            for n = m:size(EEG_list_math, 2)
                [MSC_math, f] = mscohere(EEG_list_math(:, m), EEG_list_math(:, n), hamming(L_math), round(0.5 * L_math), nfft, fs);        
                filtered_MSC_math = MSC_math(freq_indexes); % extracting the MSC values of the specific band      
                mean_MSC_math = mean(filtered_MSC_math);% computing the mean value for the band
                MSC_results_math(m,n) = mean_MSC_math; % storing the mean value in the results matrix (19x19)
                MSC_results_math(n,m) = MSC_results_math(m,n);% the results matrix is symmetrical
                
            end
        end
        MSC_tot_math(:, :, ind) = MSC_tot_math(:, :, ind) + MSC_results_math; %summing the value of MSC during arithmetics task state of this band of all subjects to then compute the mean value
        %% 

        % computing the mean value for each row for each channel to rapresent the mean coherence for each channel using topoplot function  
        % topoplot
        mean_MSC_rt = mean(MSC_results_resting, 2); 
        mean_MSC_mt = mean(MSC_results_math, 2);


        % computing tha upper part of the symmetric matrix (MSC_results_math and
        % MSC_results_resting) to represent it using the imagesc function. 
        MSC_upper_resting = MSC_results_resting;
        for p = 1:size(MSC_upper_resting, 1)
            MSC_upper_resting(p, 1:p-1) = NaN;   
        end

        MSC_upper_math = MSC_results_math;
        for q = 1:size(MSC_upper_math, 1)
            MSC_upper_math(q, 1:q-1) = NaN;  
        end
        % creating an images for each subject for each band, containing the rapresentations obtained using the imagesc function and the topoplot function, both in the resting case and in the arithmetic tasks case
        tile = tiledlayout(2,2);
        nexttile
        imagesc(MSC_upper_resting);
        colorbar;
        xticks(1:length(channels)); 
        yticks(1:length(channels)); 
        xticklabels(channels);       
        yticklabels(channels); 
        xtickangle(45);
        title('Resting');

        nexttile
        imagesc(MSC_upper_math);
        colorbar;
        xticks(1:length(channels)); 
        yticks(1:length(channels)); 
        xticklabels(channels);       
        yticklabels(channels); 
        xtickangle(45);
        title('Arithmetic task');

        nexttile
        topoplot(mean_MSC_rt, chanlocs,'colormap',jet(500));
        colorbar;
        clim([0 1]);
        title('Resting');
        
        nexttile
        topoplot(mean_MSC_mt, chanlocs,'colormap',jet(500));
        colorbar;
        clim([0 1]);
        title('Arithmetic task');

        overall_caption = append(['MSC subject ', int2str(i), ' ', band_name, ' ','band']);
        sgtitle(overall_caption);

    end    
end





%% computing and representing with topoplot and imagesc the average value of MSC of all subjects in each band
num_subjects = length(subject_indexes);
MSC_tot_resting_new = MSC_tot_resting / num_subjects;  %Dividing the sum by the number of subject so that I obtain the mean value of all subject
MSC_tot_math_new = MSC_tot_math / num_subjects;
% creating an image for each benad with the mean value of all subject,
% containing the imagesc rapresentation and the topoplot rapresentation in
% both cases (resting and arithmetics)
for ind = 1:num_bands
    figure; 
     band_name = band_names{ind}; 

     % Imagesc for resting state
    tile = tiledlayout(2,2);
    nexttile
    imagesc(MSC_tot_resting_new(:, :, ind));
    colorbar;    
    xticks(1:length(channels));
    yticks(1:length(channels));
    xticklabels(channels);
    yticklabels(channels);
    xtickangle(45);
    clim([0 1]);
    title('Resting');


    % Imagesc for arithmetic tasks
    nexttile
    imagesc(MSC_tot_math_new(:, :, ind));
    colorbar;    
    xticks(1:length(channels));
    yticks(1:length(channels));
    xticklabels(channels);
    yticklabels(channels);
    xtickangle(45);
    clim([0 1]);
    title('Arithmetic task');

    % Topoplot for resting state
    nexttile
    mean_MSC_resting = mean(MSC_tot_resting_new(:, :, ind), 2);
    mean_MSC_resting = squeeze(mean_MSC_resting); 
    topoplot(mean_MSC_resting, chanlocs, 'colormap', jet(500));    
    colorbar;
    title('Resting');
    clim([0 1]);


    % Topoplot for arithmetic task
    nexttile
    mean_MSC_math = mean(MSC_tot_math_new(:, :, ind), 2);
    mean_MSC_math = squeeze(mean_MSC_math);
    topoplot(mean_MSC_math, chanlocs, 'colormap', jet(500));
    colorbar;
    title('Arithmetic task');
    clim([0 1]);
    
    
    overall_caption = append(['Mean of all subjects band',' ', band_name]);
    sgtitle(overall_caption);

end


%% computing an image with the topoplots of the mean value of all subject in all the bands, both in resting case and in arithmetic tasks case
figure;
for ind = 1:num_bands
    band_name =band_names{ind} ;     

    % Topoplot for resting state
    subplot(3, 2, 2*ind-1); 
    mean_MSC_tot_resting = mean(MSC_tot_resting_new(:, :, ind), 2);
    topoplot(mean_MSC_tot_resting, chanlocs, 'colormap', jet(500));
    colorbar;
    clim([0 1]);
    title(['Resting band ', band_name]);

    % Topoplot for math state
    subplot(3, 2, 2*ind); 
    mean_MSC_tot_math = mean(MSC_tot_math_new(:, :, ind), 2);
    topoplot(mean_MSC_tot_math, chanlocs, 'colormap', jet(500));
    colorbar;
    clim([0 1]);
    title(['Arithmetic band ', band_name]);
end

overall_caption = append('Mean of all subjects');
sgtitle(overall_caption);




           

           
                   


function [CSA, freq, time] = computeCSA (EEG, phase)
    % parameters for the PSD using Welch method
    fs = 500;
    % defining different time windows for resting (which lasts 3 minutes)
    % and math (lasting 1 minute only)
    if phase == 1
        window_length = fs * 3;
    else
        window_length = fs * 1;
    end
    shift_ratio = 0.5;  % 50% overlap
    nfft        = 2048;
    noverlap    = fs * shift_ratio;
    % variables for CSA
    N = length(EEG);
    i_idx = 1;
    counter = 1;
    CSA = [];
    time = [];

    while i_idx + window_length < N
        % EEG windowed
        EEG_segment = detrend(EEG(i_idx:i_idx + window_length));

        %  PSD with Welch
        [PSD, freq] = pwelch(EEG_segment, hamming(fs), noverlap, nfft, fs);

        % just freq below 30
        PSD = PSD(freq < 30);

        % Save PSD in CSA matrix
        CSA(counter, :) = PSD;

        % time in the center of the window
        time(counter) = (i_idx + window_length / 2) / fs / 60;

        % move the index for next segment
        i_idx = i_idx + window_length * shift_ratio;  % shift to move the window
        counter = counter + 1;
    end
end

function time = computeTime(phase, N)
    fs = 500;
    if phase == 1
        window_length = fs * 3;
    else
        window_length = fs * 1;
    end
    i_idx = 1;
    counter = 1;
    shift_ratio = 0.5;

    while i_idx + window_length < N
        % time in the center of the window
        time(counter) = (i_idx + window_length / 2) / fs / 60;

        % move the index for next segment
        i_idx = i_idx + window_length * shift_ratio; 
        counter = counter + 1;
   end
end

