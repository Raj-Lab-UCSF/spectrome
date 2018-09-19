% list of lobe regions
% 1 = frontal, 2 = parietal, 3 = occipital, 4 = temporal, 5 = subcortical, % 6 = cerebellum
% Colors: 'k' = black, 'r' = red, 'b' = blue, 'y' = yellow, 'g' = green, 

 lobes = zeros(116,1);
 lobes(1:28) = 1;
lobes([29:42, 71:78]) = 5;
lobes(43:54) = 3;
lobes([55:56, 79:90]) = 4;
lobes(57:70) = 2;
lobes(91:116) = 6;


% new mods in limbic system
lobes([29, 30]) = 1;  % insula
lobes([31, 32]) = 1;  % ant cing
lobes([33, 34]) = 1; % mid cing
lobes([35, 36]) = 2; % post cing
lobes([37, 38]) = 4; % hippo
lobes([39, 40]) = 4; % parahippo

