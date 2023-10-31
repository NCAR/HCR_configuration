function [specLin,specDB]=getSpectra(cIQ)
% FFT and spectra
%% V
fftIQv=fft(cIQ.v,[],2);

powerRealInV=real(fftIQv).^2;
powerImagInV=imag(fftIQv).^2;

powerSignalV=powerRealInV+powerImagInV;
powerShiftedV=fftshift(powerSignalV,2);

% Reverse to get pointing direction consistent
specLin.V=fliplr(powerShiftedV);

% DB
specDB.V=10*log10(specLin.V);
%% H
if isfield(cIQ,'h')
    fftIQh=fft(cIQ.h,[],2);

    powerRealInH=real(fftIQh).^2;
    powerImagInH=imag(fftIQh).^2;

    powerSignalH=powerRealInH+powerImagInH;
    powerShiftedH=fftshift(powerSignalH,2);

    % Reverse to get pointing direction consistent
    specLin.H=fliplr(powerShiftedH);

    % DB
    specDB.H=10*log10(specLin.H);
end
end