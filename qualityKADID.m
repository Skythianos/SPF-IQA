clear all
close all

warning off

load KADID_Data2.mat

path = 'C:\Users\Public\QualityAssessment\KADID-10k\images';

numberOfImages = size(dmos, 1);
Scores = zeros(numberOfImages, 1);

Features = zeros(numberOfImages, 57);

parfor i=1:numberOfImages
    if(mod(i,1000)==0)
        disp(i);
    end
    img  = imread( char(strcat(path, filesep, string(dist_img(i)))) );
    Features(i,:) = getFeatures(img);
end

PLCC = zeros(1,20);
SROCC= zeros(1,20);

for i=1:20
    disp(i);
    [Train, Test] = splitTrainTest(dist_img);

    TrainFeatures = Features(Train,:);
    TestFeatures  = Features(Test,:);
    
    YTest = dmos(Test);
    YTrain= dmos(Train);

    Mdl = fitrgp(TrainFeatures,YTrain,'KernelFunction','rationalquadratic','Standardize',true);
    Pred= predict(Mdl,TestFeatures);
    
    PLCC(i) = round(corr(Pred, YTest),3);
    SROCC(i)= round(corr(Pred, YTest, 'Type', 'Spearman'),3);
end

disp('----------------------------------');
X = ['Average PLCC after 20 random train-test splits: ', num2str(round(mean(PLCC(:)),3))];
disp(X);
X = ['Average SROCC after 20 random train-test splits: ', num2str(round(mean(SROCC(:)),3))];
disp(X);