function [PC] = getPhaseCongruencyImage(I)

    [Y,X,ch] = size(I);

    if(ch==3)
        I = rgb2gray(I);
    end
    
    % This time we have to use exactly two scales, as this is required by
    % Felsberg's phase congruency method. As he suggests, let's use
    % three-ocatave filters and leave a three-octave spacing between them.
    % We want small filters to pick out the fine detail in this image.
    cw = [3,24];

    % Construct new filters, as before
    filtStruct = createMonogenicFilters(Y,X,cw,'lg',0.41);

    % Find monogenic signal, as before
    [m1,m2,m3] = monogenicSignal(I,filtStruct);

    % Now use the phase congruency algorithm. The fourth parameter is a
    % threshold between 0 and 1 used for noise supression. You will always need
    % to use this to get reasonable results. Somewhere between 0 and 0.1 should
    % do in most cases.
    PC = phaseCongruency(m1,m2,m3,0.05);

end

