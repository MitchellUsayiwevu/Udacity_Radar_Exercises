close all;
clear all;
clc;

%% Radar Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = 77*10^9;
d_max = 200;
d_res = 1;
v_max = 100;
c = 3e8;

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant

target_pos = 50;
target_vel = 23;

%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
B = c/(2*d_res);
% chirp using the requirements above.
T_chirp = (5.5*2*d_max)/c;
S = B/T_chirp;
%Operating carrier frequency of Radar
fc= 77e9;             %carrier freq


%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation.
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp.
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*T_chirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
t_d=zeros(1,length(t));
delta_t = t(2)-t(1);


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.

for i=1:length(t)


    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity.
   if i>1
     r_t(i) =  r_t(i-1) + (delta_t*target_vel);
   else
      r_t(i) =  target_pos + (delta_t*target_vel);
   end
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal.
    t_d(i) =  (2*r_t(i))/c;
    Tx(i) = cos(2*pi*((fc*i)+((S*i^2)/2)));
    Rx (i)  = cos(2*pi*((fc*(i-t_d(i)))+((S*(i-t_d(i))^2)/2)));

    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);

end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix1d = reshape(Mix, [Nr, Nd]);
 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.

signal_fft = fft(Mix1d,Nr,1);
 % *%TODO* :
% Take the absolute value of FFT output
signal_fft_abs = abs(signal_fft/Nr);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
P1 = signal_fft_abs(1:Nr/2 +1);

%plotting the range
 % *%TODO* :
 % plot FFT output

figure()
plot(0:(Nr/2), P1)
xlabel('Range bin')
ylabel('|Amplitude|')
title('Range FFT Output')
grid on

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.

Tr = 8;
Td = 4;
Gr = 4;
Gd = 2;
CUT = 1;
offset = 13;


% decent
%Tr = 12;
%Td = 12;
%Gr = 4;
%Gd = 4;
%CUT = 1;
%offset = 10;



% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation

% *%TODO* :
% offset the threshold by SNR value in dB


% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);
threshhold_val = zeros(1,1);

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

Threshhold_block = zeros(Nr/2,Nd);
training_cells_total =  (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1);

for i = (Tr+Gr+1) : ((Nr/2) - (Tr+Gr+1))
  for j = (Td+Gd+1) : (Nd - (Td+Gd+1))

   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR


     % loop to get to each of the training cells.
     % indeces are relative to the CUT index, i and j
    noise_level_sum = 0;
     for k = (i -(Tr+Gr)) : (i + (Tr+Gr))

       for l = (j - (Td+Gd)) : (j + (Td+Gd))

         if k >= (i-Gr) && k <= (i + Gr) && l >=(j-Gd) && l <= (j + Gd)
           continue;
         else
           noise_level = RDM(k,l);
%           noise_level_pow = db2pow(noise_level);
           noise_level_pow = 10.^(noise_level / 10);
           noise_level_sum = noise_level_sum + noise_level_pow;
         end

       end
     end

    noise_level_avg = noise_level_sum/ training_cells_total;
%     noise_level_avg_db =  pow2db(noise_level_avg);
     noise_level_avg_db = 10 * log10(noise_level_avg);
     threshhold_val = noise_level_avg_db + offset;

     if RDM(i,j)>threshhold_val
       Threshhold_block(i,j) = 1;
    else
      Threshhold_block(i,j) = 0;
    end


  end
end




% *%TODO* :
% The process above will generate a thresholded block, which is smaller
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0.



% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.


figure;
surf(doppler_axis,range_axis, Threshhold_block);
colorbar;

%%%

##This section uses code from mentor who helped me debug my code using convolution
## Helped with trying out different parameter combinations quickly as this has a faster execution time compared to the code I wrote with nested for loops

%pow2db = @(x) 10*log10(x);
%db2pow = @(x) 10.^(x / 10);
%mask = ones(2 * Tr + 2 * Gr + 1, 2 * Td + 2 * Gd + 1);
%centre_coord = [Tr + Gr + 1, Td + Gd + 1];
%mask(centre_coord(1) - Gr : centre_coord(1) + Gr, centre_coord(2) - Gd : centre_coord(2) + Gd) = 0;
%mask = mask / sum(mask(:));
%% Convolve, then convert back to dB to add the offset
%% The convolution defines the threshold
%threshold = conv2(db2pow(RDM), mask, 'same');
%threshold = pow2db(threshold) + offset;
%% Any values less than the threshold are 0, else 1
%RDM(RDM < threshold) = 0;
%RDM(RDM >= threshold) = 1;
%
%figure;
%surf(doppler_axis,range_axis, RDM);
%colorbar;

%%



