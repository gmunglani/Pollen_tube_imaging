function threshold = distribution_bimodal(BF1,thresh_sd)

x = [];
x = reshape(BF1,1,size(BF1,1)*size(BF1,2));
x = double(x);

%pd = fitdist(x','Normal');
%threshold1 = round(pd.mu+3*pd.sigma)

pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
                         p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
                     
pStart = .5;
muStart = quantile(x,[.1 .90]);
sigmaStart = sqrt(var(x) - .25*diff(muStart).^2);
start = [pStart muStart sigmaStart sigmaStart];

lb = [0 -Inf -Inf 0 0];
ub = [1 Inf Inf Inf Inf];

statset('mlecustom');

options = statset('MaxIter',300, 'MaxFunEvals',600);
paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start, ...
                          'lower',lb, 'upper',ub, 'options',options)

threshold = round(real(paramEsts(2))+real(thresh_sd*paramEsts(4)))

% figure
% bins = 0:10:max(x);
% h = bar(bins,histc(x,bins)/(length(x)*.5),'histc');
% h.FaceColor = [.9 .9 .9];
% xgrid = linspace(1.1*min(x),1.1*max(x),200);
% pdfgrid = pdf_normmixture(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
% hold on
% plot(xgrid,pdfgrid,'-')
% hold off
% xlabel('x')
% ylabel('Probability Density')                  