function [colorFeatures] = getColorStatistics(imgRGB)

    if(isa(imgRGB,'uint8'))
        imgRGB = im2double(imgRGB);
    end

    R = imgRGB(:,:,1);
    G = imgRGB(:,:,2);
    B = imgRGB(:,:,3);
    
    [M, N, ~] = size(imgRGB);
    
    for i=1:M
        for j=1:N
            if(R(i,j)==0)
                R(i,j)=0.0001;
            end
            if(G(i,j)==0)
                G(i,j)=0.0001;
            end
            if(B(i,j)==0)
                B(i,j)=0.0001;
            end
        end
    end
    
    logR = log(R);
    logG = log(G);
    logB = log(B);
    
    R1 = logR - mean(logR(:));
    G1 = logG - mean(logG(:));
    B1 = logB - mean(logB(:));
    
    l1 = (R1 + G1+ B1)/(sqrt(3));
    l2 = (R1 + G1 - 2*B1)/(sqrt(6));
    l3 = (R1 - G1)/(sqrt(2));
    
    pd1 = fitdist(l1(:),'Normal');
    pd2 = fitdist(l2(:),'Normal');
    pd3 = fitdist(l3(:),'Normal');
    
    colorFeatures = [pd1.mu, pd1.sigma, pd2.mu, pd2.sigma, pd3.mu, pd3.sigma];

end

