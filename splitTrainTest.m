function [Train,Test] = splitTrainTest(Names)

    numberOfImages = size(Names,1);
   
    Train = false(numberOfImages,1);
    Test  = false(numberOfImages,1);
    
    p = randperm(81);
    
    train = p(1:65); 
    
    for i=1:numberOfImages
        name = Names{i};
        tmp = char(name);
        tmp = str2double(tmp(2:3));
        if( ismember(tmp,train) )
            Train(i)=true;
        else
            Test(i)=true;
        end
    end

end
