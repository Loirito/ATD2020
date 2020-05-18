function STFT(fileData)
    %fs = 799.1514;
    fs = 50;
    Ts = 1/fs;
    data = fileData{1};
    N = numel(data);
    t = N*Ts;
    Tframe = 0.005*t;
    
    Toverlap = Tframe/2;
    
    Nframe = round(Tframe*fs);
    h = hamming(Nframe);
    
    Noverlap = round(Toverlap*fs);
    freq_relev = [];
    f = linspace(-fs/2,fs/2,Nframe);
    x =  find(f>=0);
    for ii = 1:Nframe-Noverlap:N-Nframe
        x_frame = data(ii:ii+Nframe-1).*h;
        m_X_frame = abs(fftshift(fft(x_frame)));
        freq_relev = horzcat(freq_relev,m_X_frame(x));
    
    end
    figure()
    spectrogram(data,Nframe,Noverlap,[],fs,'yaxis')

    [s,f,t,p] = spectrogram(data,Nframe,Noverlap,[],fs);
    
end