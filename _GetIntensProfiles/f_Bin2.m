%% The task of the program is to do software binning 2 on images
function [ResImage] = f_Bin2(InitImage)
s = size(InitImage);
ResImage = zeros(s/2);
for Lin = 1:2:s(1) 
    for Col = 1:2:s(2)
        ResImage(floor(Lin/2) + 1, floor(Col/2) + 1) = sum(InitImage(Lin:Lin + 1, Col)...
            + InitImage(Lin:Lin + 1, Col + 1));
    end
end
ResImage(find(ResImage)) = 1;
figure, imshow(InitImage, []);
figure, imshow(ResImage, []);
