function [MSCN_img]= MSCN(img)

    if(size(img,3)==3)
        img = rgb2gray(img);
    end

    if(isa(img,'uint8'))
        img = im2double(img);
    end

    window = fspecial('gaussian',7,7/6); 
    window = window/sum(sum(window));
   
    mu = filter2(window, img, 'same');
    mu_sq = mu.*mu;

    sigma = sqrt(abs(filter2(window, img.*img, 'same') - mu_sq));
    imgMinusMu = (img-mu);

    MSCN_img =imgMinusMu./(sigma +1);
end