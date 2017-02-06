%% Starting codes for the part-III in the HA2 of CS596, Fall, 2016.

% Fill in the codes between "%PLACEHOLDER#start" and "PLACEHOLDER#end"

%% step 1: generate dataset that includes both positive and negative samples,
% where each sample is described with two features. 
%250 samples in total.

file_x = 'adhd_x.xlsx';
file_y = 'adhd_y.xlsx';

X_dataM = xlsread(file_x,'A1:ABB1065');
Y_labels = xlsread(file_y,'A1:ABB1');
Y_labels = Y_labels';
Y_labels(Y_labels > 0) = 1;

nsamples  = 700; 
ntrain    = nsamples*4/5;  
nfeatures = [10,100,200,300,400,500,600,700,800,900,1000];
fidx      = 1;
avg_error   = zeros(fidx,1);
std_error   = zeros(fidx,1);

%[X,y]=getDataset(); % note that y contains only 1s and 0s,

   
for i=1:fidx
    
   %plotting all samples
   idx1 = find(Y_labels == 0); % object indices for the 1st class
   idx2 = find(Y_labels > 0);
   h=subplot(1,3,1);  hold on;
    
    %X_data = X_dataM(1:nfeatures(i),:);
    X_data = X_dataM(6,:);
    size(X_data)
    
    % no more variables are needed
    %plot(X_data(idx1,1), X_data(idx1,2), 'r*');
    plot(X_data(1:nfeatures(i),idx1), X_data(1:nfeatures(i),idx1), 'r*');
    %plot(X_data(idx2,1), X_data(idx2,2), 'b*');
    plot(X_data(1:nfeatures(i),idx2), X_data(1:nfeatures(i),idx2), 'b*');
    %axis tight
    xlabel('x_1');
    ylabel('x_2');
    title('All samples');
    
%number of training samples



 %PLACEHOLDER#start
 % write you own code to randomly pick up 100 samples for training and use
 % the rest for testing. 
 %idx = round((249)*rand(249,1));
 
 
 
 randz     = randperm(nsamples);
 rand1     = randz(1:ntrain);
 rand2     = randz(ntrain+1:end);


 %trainX=X_data(rand1, :); % training samples,
 trainX=X_data(:,rand1); % training samples,
 trainy=Y_labels(rand1); % labels of training samples    
 testX=X_data(:,rand2); % testing samples
 testy=Y_labels(rand2); % labels of testing samples   
    
 %PLACEHOLDER#end
 
 % plot the samples you have pickup for training, check to confirm that both negative
 % and positive samples are included. 
 
    idx1 = find(trainy == 0); % object indices for the 1st class
    idx2 = find(trainy > 0);
    h=subplot(1,3,2); hold on;
    % no more variables are needed
    
    %plot(trainX(idx1,1), trainX(idx1,2), 'r*');
    %plot(trainX(idx2,1), trainX(idx2,2), 'b*');
    
    plot(trainX(1:nfeatures(i),idx1), trainX(1:nfeatures(i),idx1), 'r*');
    plot(trainX(1:nfeatures(i),idx2), trainX(1:nfeatures(i),idx2), 'b*');
    
    %axis tight
    
    xlabel('x_1');
    ylabel('x_2');
    title('training samples');

    idx1 = find(testy == 0); % object indices for the 1st class
    idx2 = find(testy == 1);
    % no more variables are needed
    h=subplot(1,3,3); hold on; 
    
    %plot(testX(idx1,1), testX(idx1,2), 'r*');
    %plot(testX(idx2,1), testX(idx2,2), 'b*');
    
    plot(testX(1:nfeatures(i),idx1), testX(1:nfeatures(i),idx1), 'r*');
    plot(testX(1:nfeatures(i),idx2), testX(1:nfeatures(i),idx2), 'b*');
    
    %axis tight
    xlabel('x_1');
    ylabel('x_2');
    title('testing samples');    
    
%% step 2: train logistic regression model with the gradient descent method 

bHat = glmfit(trainX',trainy,'binomial');

%% step 3: Use the model to get class labels of testing samples.
 %PLACEHOLDER#start
    % with the learned model, apply the logistic model over testing
    % samples; hatProb is the probability of belonging to the class 1.
    %hatProb = rand(length(testy),1);  % variant of classification
    testX = testX';
    hatProb = 1./(1+exp( -[ones( size(testX,1),1 ), testX] *bHat));

 %PLACEHOLDER#end
 
 % predict the class labels with a threshold
haty=(hatProb>=0.5); 
%% step 4: evaluation
 
   %compare haty and testy to calculate average error and standard deviation  %   
   avgErr       = mean(abs(haty-testy));
   stdErr       = std(abs(haty-testy));
   avg_error(i) = mean(abs(haty-testy));
   std_error(i) = std(abs(haty-testy));
 
end  
 
 fprintf('average error:%f (%f)\n', avgErr, stdErr);
 %plot(nfeatures,t_error);
 