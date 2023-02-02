function [Thresh,error_subspace]=Calculate_Alpha(Data,labels,Vecmode,Curveoption)
%% This piece of code runs the PLDA function with different 
%  values of alpha and plot two error curves
%  The projection metric distance of two consequent subspaces
%  For    
%% Load the data set
[Nfeatures,Npnt] = size(Data);
%%
nPLDA=Curveoption.nPLDA; 
counter=0;
x=Curveoption.low:Curveoption.step:Curveoption.high;
for Alpha=x
    Alpha;
    counter=counter+1; 
    [Vec(:,:,counter)] = PLDA_in_PCAspace(Data,labels,Alpha,nPLDA,Vecmode);
    if counter>1
          Vec(:,:,counter)=Vec(:,:,counter)*diag(sign(diag(Vec(:,:,counter)'*Vec(:,:,counter-1))));               
          error_subspace(counter-1)= Projection_metric( Vec(:,:,counter),Vec(:,:,counter-1)); 
    end
end 
%% Calculate twice the half life of Alpha,

%Fit an exponential to log(error_subspace)
f =@(a,b,x) log(a)-b*x; 
options = fitoptions('Method','linearLeastSquares');
F_fitted = fit(x(1:end-1)',log(error_subspace)', f, ...
    'StartPoint', [1,x(1)], ...
    'Lower', [0,0],'Robust','LAR');
coeff=coeffvalues(F_fitted);%Get coefficients of fitted function
Thresh=log(2)*(1/coeff(2));%Calculate twice the half life = 2(log(2)/b)

% figure
% plot(x(1:end-1),error_subspace,'linewidth',2)
% title('Stability of subspace w.r.t. Alpha','fontsize',24)
% ylabel('Projection metric between two consequent subspaces','fontsize',20)
% xlabel('Alpha','fontsize',20)
% set(gca,'fontsize',20)
% yL = get(gca,'YLim');
% hold on
% line([Thresh Thresh],yL,'Color','r','linewidth',2);
% grid on
