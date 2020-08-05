close all;
% p_ResFit = 0.98;

AverageIntensBright = CytMaxima1552;

NbBinsAvIntens = ceil(max(AverageIntensBright) - min(AverageIntensBright)) / 2.5;
figure(1),  
[Nb,XOut] = hist(AverageIntensBright(:, 1), NbBinsAvIntens);
bar(XOut, Nb);
% title('DB1665,gain70,exp999ms,bin1,13planes,05space,date:070503');      
% xlabel('MT maxima');
% ylabel('N');
grid on

% x_Fit = XOut(1):0.05:XOut(length(XOut));
% Spl = csaps(XOut, Nb, p_ResFit, x_Fit);
% 
% figure(2);
% hold on;
% bar(XOut, Nb);
% %ylabel('N');
% %figure(2);
% % plot(NbsHist, 'bs', 'MarkerSize', 3);
% % hold on
% plot(x_Fit, Spl, '-r', 'LineWidth', 2.5);
% hold off
% grid on
% title('DB1665,gain70,exp999ms,bin1,13planes,05space,date:070503');      
% xlabel('MT maxima');
% ylabel('N');
% 
% [DistrMax, index] = simple_max(Spl);
% DistrMaxPos = x_Fit(index);
% Std = std(Spl);
% %HalfHeightWidth = Spl(); 
% 
% save(FileName);

% 
% 
% NbBins = ceil(simple_max(Result));
% [NbsHist, Xout] = hist(Result, NbBins);
% 
% x_Fit = 1:0.2:length(NbsHist);
% Spl = csaps(Xout, NbsHist, p_ResFit, Xout);
% 
% figure(1);
% hold on;
% bar(Xout, NbsHist);
% %ylabel('N');
% 
% %figure(2);
% % plot(NbsHist, 'bs', 'MarkerSize', 3);
% % hold on
% plot(Xout, Spl, '-r', 'LineWidth', 2.5);
% hold off
% grid on
% ylabel('N');