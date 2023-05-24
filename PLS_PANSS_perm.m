clear;clc
load('E:\NAS\LapTop_ssj\SZ-Gen\Writing\Submit\BMCmedicine\Revise1\Response2\GitHub\demo\input\Ydata.mat') % load dependent vadiables
% z-score:
Y=zscore(Y);
load('E:\NAS\LapTop_ssj\SZ-Gen\Writing\Submit\BMCmedicine\Revise1\Response2\GitHub\demo\input\Xdata.mat') % load independent variables
X=parcelExpression; X(isnan(X))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%perform full PLS and plot variance in Y explained by top 15 components
%typically top 2 or 3 components will explain a large part of the variance
%(hopefully!)
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
%save('MCM_fun_pls_result', 'XL', 'YL', 'XS', 'YS', 'BETA', 'PCTVAR', 'MSE', 'stats');
dim=15;
plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
grid on
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim); % no need to do this but it keeps outputs tidy
%%% plot correlation of PLS component 1 with dependent variables:
figure
plot(XS(:,1),Y,'r.')
[R,p]=corrcoef(XS(:,1), Y);
xlabel('XS scores for PLS component 1','FontSize',14);
ylabel('Cobre t-statistic- lh','FontSize',14);
grid on

rep=1000;
p_perm=[]; Rsq_perm=[]; 
for dim=1:8
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    temp=cumsum(100*PCTVAR(2,1:dim));
    Rsquared = temp(dim);
    for j=1:rep
        %j
        order=randperm(size(Y,1));
        Yp=Y(order,:);
        
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);
        
        temp=cumsum(100*PCTVAR(2,1:dim));
        Rsq_perm(j) = temp(dim);
    end
    p_perm(dim)=length(find(Rsq_perm>=Rsquared))/rep;
end
figure
plot(1:dim, p_perm,'ok','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Number of PLS components','FontSize',14);
ylabel('p-value','FontSize',14);
grid on