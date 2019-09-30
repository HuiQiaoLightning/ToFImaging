% ***Show reconstructed transient image***
%
% see our paper: F. Heide, M. Hullin, J. Gregson, W. Heidrich.: 
% "Low-budget Transient Imaging using Photonic Mixer Devices.", SIGGRAPH 2013.
%

% Copyright (C) 2013. Felix Heide
% Email: fheide@cs.ubc.ca

%Clear workspace
clear all
close all hidden

%Dataset
I_reconstructed = load('./reconstruction_I.mat');

%Scale and gamma
scaling=1.0;

I_rec = I_reconstructed.I_reconstructed;
I_rec = max(I_rec,0);
I_rec = I_rec./max(I_rec(:));
I_rec = scaling * max(min(I_rec,1),0);

I_rec = I_rec .^ (1/2.2); %Gamma
scalefactor = 2;

%Show transient image and reconstruction
showrecfig = figure();
for t = 1:1:size(I_rec,3) - 100

    I_show = imresize( I_rec(:,:,t), scalefactor,'bilinear');
    
    %Flip
    I_show = I_show(end:-1:1, :);
    
    %Show image
    imshow( I_show ), title(sprintf('Rec Frame %d', t));
    pause(0.03)
    
end

fprintf('\nDone\n')
