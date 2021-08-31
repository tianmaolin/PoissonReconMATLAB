function imO = wihteBoundary(im)
im(isnan(im)) = 1;
label = bwlabel(im);
imO = im;
imO(label == 1) = 0;