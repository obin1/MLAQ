clear all; close all;
% Load data
JC = readtable('C:\Users\islam\OneDrive - Johns Hopkins\Wexler\Atmos\FormO3NOx-master\C_and_J.txt');
S = readtable('C:\Users\islam\OneDrive - Johns Hopkins\Wexler\Atmos\FormO3NOx-master\S.txt');
JC(JC.Time_min_==1440,:)=[]; %remove last concentrations
JC = JC(:,3:14);
S = S(:,3:12);

A = [0,  1, -1,  0,  0,  0,  0,  0,  0,  0;  % O3
     1,  0, -1,  0,  0,  0, -1,  0,  0,  0;  % NO
     -1,  0,  1,  0,  0,  0,  1, -1,  0,  0;  % NO2
     0,  0,  0, -1, -1, -1,  0,  0,  0,  0; % HCHO
     0,  0,  0,  2,  0,  1, -1,  0,  0,  1; % HO2.
     0,  0,  0,  0,  0,  0,  0,  0, -1, -1;];  % HO2H
 
% Visualize
scatter(table2array(JC(:,1)), table2array(S(:,1)),'.'); title('R1');...   
scatter(table2array(JC(:,3)), table2array(S(:,1)),'.');
scatter(table2array(JC(:,4)), table2array(S(:,1)),'.');

scatter3(table2array(JC(:,3)),table2array(JC(:,4)),table2array(S(:,1)),'.');...
    title('R1'); xlabel('NO'); ylabel('NO2'); zlabel('S');  
scatter3(table2array(JC(:,1)),table2array(JC(:,2)),table2array(S(:,2)),'.');...
    title('R2'); xlabel('Light'); ylabel('O3'); zlabel('S');  

% Test/Train Split
[trainInd,valInd,testInd] = dividerand(size(JC, 1),.9,0,.1);
JC_train = JC(trainInd,:);
JC_test = JC(testInd,:);
S_train = S(trainInd,:);
S_test = S(testInd,:);

% Store relevant species for each reaction
species = {};
for r = 1:10
    species{r} = [1; find(A(:,r))+1]; 
end

% Train SVM Models
svm = {};
for r = 1:10
    svm{r} = fitrsvm(JC_train(:,species{r}), S_train(:,r),'KernelFunction','gaussian');
end
save('svmr.mat');

abs_err = zeros(10,1); rel_err = zeros(10,1);
% Test predictions
for r = 1:10
    preds{r} = predict(svm{r}, JC_test(:,species{r}));
    abs_err(r) = immse(preds{r}, table2array(S_test(:,r)));
    rel_err(r) = abs_err(r) / mean(table2array(S_test(:,r)));
end

% Plot results
scatter(1:10, abs_err); xlabel('R'); ylabel('MSE');
scatter(1:10, rel_err * 100); xlabel('R'); ylabel('Rel. MSE (%)');



