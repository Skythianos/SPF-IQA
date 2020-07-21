function [featureVector] = getFeatures(imgRGB)

imgGray = rgb2gray(imgRGB);

% 
FD = getFD(imgRGB);
[featureVector_1,~] = histcounts(FD(:), [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3], 'Normalization', 'probability');

% 
[~,cH,cV,cD]=dwt2(imgRGB,'sym4','mode','per');
% 
featureVector_2 = getFirstDigitDistribution(cH);

featureVector_3 = getFirstDigitDistribution(cV);

featureVector_4 = getFirstDigitDistribution(cD);
% 
[Mag,~] = imgradient(imgGray);
featureVector_5 = getFirstDigitDistribution(Mag);
% 
PC = getPhaseCongruencyImage(imgRGB);
% 
colorFeatures = getColorStatistics(imgRGB);
% 
featureVector = [featureVector_1, featureVector_2, featureVector_3, featureVector_4, featureVector_5, colorFeatures,...
    getColorfulness(imgRGB), getGlobalContrastFactor(imgRGB),getDarkChannelFeature(imgRGB), entropy(imgRGB), mean(PC(:)), ...
    ];
%,...
    %entropy(PC), entropy(imgRGB), getGlobalContrastFactor(imgRGB),...
    %getColorfulness(imgRGB), getDarkChannelFeature(imgRGB), getColorStatistics(imgRGB)];

%     [w, h, ch] = size(img);
%     if(ch==3)
%         img=rgb2gray(img);
%     end
%         
%     GCD = gcd(w,h);
%     
%     w_width = GCD;
%     w_height = GCD;
%     stepSize = GCD/2;
%     
%     featureVector = [];
%     
%     for i=1:stepSize:w-w_width
%         for j=1:stepSize:h-w_height
%             window = img(i:i+w_width, j:j+w_height);
%             
%             FD = getFD(window);
%             [N,~]=histcounts(FD(:), [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3], 'Normalization', 'probability');
%             
%             featureVector = [featureVector, N];
%         end
%     end
    
end
