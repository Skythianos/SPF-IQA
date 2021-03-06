function LE = localEnergy(m1,m2,m3,wl_ind)
%
%   LE = localEnergy(m1,m2,m3,wl_ind)
%
% Calculates the local energy from a monogenic signal (m1,m2,m3) of an
% image or video. The phase is calculated for each scale independently.
%
% Alternatively, wl_ind, a parameter selecting the index of the wavelength
% to use may be passed.
%
% Chris Bridge, Institute of Biomedical Engineering, University of Oxford
% christopher.bridge@eng.ox.ac.uk

if nargin < 4
   LE = m1.^2 + m2.^2 + m3.^2;
else
   LE = m1(:,:,:,wl_ind).^2 + m2(:,:,:,wl_ind).^2 + m3(:,:,:,wl_ind).^2;
end