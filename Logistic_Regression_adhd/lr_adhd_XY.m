
file_x = 'adhd_x.xlsx';
file_y = 'adhd_y_C.xlsx';

X_dataM = xlsread(file_x, 'A1:ABB500');
Y_labels = xlsread(file_y,'A1:ABB1');
Y_labels = Y_labels';
Y_labels(Y_labels > 0) = 1;


ntrain      = 10;  
nfeatures   = 500;

idx1 = find(Y_labels==0);
idx2 = find(Y_labels>0);
h=subplot(1,3,1);  hold on;
X_data  = X_dataM(1:nfeatures,:);
size(X_data)
    
plot(X_data(:,idx1), X_data(:,idx1), 'r*');
plot(X_data(:,idx2), X_data(:,idx2), 'g*');    
xlabel('x_1');
ylabel('x_2');
title('All samples');

 nsamples    = 730; 
 rand   = randperm(nsamples);
 r1     = rand(1:ntrain);
 r2     = rand(ntrain+1:end);
 trainX=X_data(:,r1); 
 trainy=Y_labels(r1); 
 testX=X_data(:,r2); 
 testy=Y_labels(r2); 
 
idx1  = find(not(trainy));
idx2 = find(trainy);
h=subplot(1,3,2); hold on; 
plot(trainX(:,idx1), trainX(:,idx1), 'r*');
plot(trainX(:,idx2), trainX(:,idx2), 'g*');

xlabel('x_1');
ylabel('x_2');
title('training samples');
idx1  = find(not(testy));
idx2  = find(testy);

h=subplot(1,3,3); hold on; 
    
plot(testX(:,idx1), testX(:,idx1), 'r*');
plot(testX(:,idx2), testX(:,idx2), 'g*');
xlabel('x_1');
ylabel('x_2');
title('testing samples');    

bHat = glmfit(trainX',trainy,'binomial');
testX = testX';
hatProb = 1./(1+exp( -[ones( size(testX,1),1 ), testX] *bHat));
haty=(hatProb>=0.5); 
 
avgErr       = mean(abs(haty-testy));
stdErr       = std(abs(haty-testy));
   
fprintf('average error:%f, std error:%f\n', avgErr, stdErr);
 
 