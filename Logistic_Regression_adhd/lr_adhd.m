%% CS596:   Final Project
%% Author:  Kim Beswick
%% Project: ADHD Dataset Evaluation using Logistical Regression 

%file_x = 'adhd_clean.xlsx';
file_x = 'adhd_clean_w1st4.xlsx';
file_y = 'adhd_y.xlsx';

X_dataM = xlsread(file_x, 'A1:ABB1065');
%X_dataM  = xlsread(file_x, 'A1:ABB4');
Y_labels = xlsread(file_y,'A1:ABB1');
Y_labels = Y_labels';
Y_labels(Y_labels > 0) = 1;

nsamples       = 730; 

n              = 3;
range          = 1:n;
best50features = zeros(1,50);
ratio          = [1/2, 2/3, 3/4, 4/5];
%nfeatures   = 10*range;
nfeatures   = 100;
%nfeatures   = 1:n;
D =            {'normal','binomial','poisson'};
%str            = string(D);               

fidx        = n; % number of iterations to increase number of features
avg_error   = zeros(fidx,1);
std_error   = zeros(fidx,1);
start       = 1;
min_avg     = 9999;
min_std     = 9999;
minidx      = 9999;


for i=1:fidx
   %crt = ratio(i);
   %ntrain       = nsamples*crt;  
   
   %plotting all samples
   %nfeatures   = randperm(1057,50);
 
   
   %idx1 = find(Y_labels == 0); % object indices for the 1st class
   idx1a = find(not(Y_labels));
   
   %idx2 = find(Y_labels == 1);  % second class
   idx2a = find(Y_labels);
   
   h=subplot(1,3,1);  hold on;
    
    %X_data = X_dataM(1:nfeatures(i),:);
    
    % X data range 
    X_data  = X_dataM(1:nfeatures,:);
    %X_data   = X_dataM((nfeatures(i)-9):nfeatures(i),:);
    
    %X_data(X_data > 0) = 1;
    
    %f_subj         = sum(X_data(1,:));
    %m_subj         = sum(X_data(2,:));
    
    %female_overlap = and(X_data(1,:), Y_labels'); 
    %female_adhd     = sum(female_overlap);
    %male_overlap   = and(X_data(2,:), Y_labels'); 
    %male_adhd       = sum(male_overlap);

    %f_percent      = female_adhd/f_subj;
    %m_percent      = male_adhd/m_subj;
    
    size(X_data)
    
    % no more variables are needed
    %plot(X_data(1:nfeatures(i),idx1), X_data(1:nfeatures(i),idx1), 'r*');
    plot(X_data(:,idx1a), X_data(:,idx1a), 'b+');
    %plot(m_percent, m_percent, 'r*');
    
    %plot(X_data(1:nfeatures(i),idx2), X_data(1:nfeatures(i),idx2), 'b*');
    plot(X_data(:,idx2a), X_data(:,idx2a), 'g+');      %'Color',[.9,.5,.4]);
    %plot(f_percent, f_percent, 'b*');
    
    %axis tight
    xlabel('x_1');
    ylabel('y');
    title('All samples');
    
 randz     = randperm(nsamples);
 rand1     = randz(1:ntrain);
 rand2     = randz(ntrain+1:end);


 %trainX=X_data(rand1, :); % training samples,
 trainX=X_data(:,rand1); % training samples,
 trainy=Y_labels(rand1); % labels of training samples    
 testX=X_data(:,rand2); % testing samples
 testy=Y_labels(rand2); % labels of testing samples   

    %idx1 = find(trainy == 0); % object indices for the 1st class
    idx1b  = find(not(trainy));
    %idx2 = find(trainy > 0);
    idx2b = find(trainy);
    
    h=subplot(1,3,2); hold on;
    % no more variables are needed
    
    %plot(trainX(1:nfeatures(i),idx1), trainX(1:nfeatures(i),idx1), 'r*');
    %plot(trainX(1:nfeatures(i),idx2), trainX(1:nfeatures(i),idx2), 'b*');
    
    plot(trainX(:,idx1b), trainX(:,idx1b), 'b+');
    plot(trainX(:,idx2b), trainX(:,idx2b), 'g+');
    
    %axis tight
    
    xlabel('x_1');
    ylabel('y');
    title('training samples');
    
    idx1c  = find(not(testy));
    idx2c  = find(testy);
    
    % no more variables are needed
    h=subplot(1,3,3); hold on; 
    
    %plot(testX(idx1,1), testX(idx1,2), 'r*');
    %plot(testX(idx2,1), testX(idx2,2), 'b*');
    
    %plot(testX(1:nfeatures(i),idx1), testX(1:nfeatures(i),idx1), 'r*');
    %plot(testX(1:nfeatures(i),idx2), testX(1:nfeatures(i),idx2), 'b*');
    
    plot(testX(:,idx1c), testX(:,idx1c), 'b+');
    plot(testX(:,idx2c), testX(:,idx2c), 'g+');
    
    xlabel('x_1');
    ylabel('y');
    title('testing samples');    

   bHat = glmfit(trainX',trainy,D{i});  %,'link', 'probit');

   testX = testX';
   hatProb = 1./(1+exp( -[ones( size(testX,1),1 ), testX] *bHat));

   haty=(hatProb>=0.5); 
  
   avgErr       = mean(abs(haty-testy));
   stdErr       = std(abs(haty-testy));
   
   
   avg_error(i) = mean(abs(haty-testy));
   std_error(i) = std(abs(haty-testy));
   if (avg_error(i) < min_avg)
       min_avg = avg_error(i);
       min_std = std_error(i);
       minidx = i;
       best50features = nfeatures;
   end
   
 
end  
 
 fprintf('average error:%f (%f)\n', avgErr, stdErr);
 %plot(nfeatures,t_error);
 %plot(nfeatures(1:fidx), avg_error);
 %x = cellstr(D);
 %x = strings(1,3);
 %v(1) = D{1};
 %v(2) = D{2};
 %v(3) = D{3};
 x = 1:3;
 
 %plot(std_error);
 %set(gca,'xticklabel',x')
 %plot(nfeatures(1:fidx), std_error);
 plot(x, std_error);
 hold on;
 plot(x, avg_error);
 %plot(D{}, avg_error);
 %plot(nfeatures(1:fidx), avg_error);
 %plot(ratio, avg_error);
 xlabel('Normal                              Binomial                         Poisson');
 ylabel('Error');
 title('Normal, Binomial and Poisson Distribution Error');    
 legend('Standard Error','Average Error');
 
 %h = histfit(avg_error);
 %hold on;
 %histfit(std_error);
 
 