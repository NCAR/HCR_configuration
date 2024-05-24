% Test the limits of the standard aircraft broadening width correction

clear all
close all

delta=0.1:0.1:1;

widthRaw=0:0.1:5;

figure('Renderer','painters')
hold on

for ii=1:length(delta)
    squares=widthRaw.^2-delta(ii).^2;
    squares(squares<0)=0.1;

    widthCorr=sqrt(squares);

    plot(widthCorr,'-','LineWidth',1.5)
end

grid on
box on