function [Z_to_test_angle_im, z_vct_full] = create_3D_saw_tooth(base_half_width, angles)

angles_of_base = (180 - angles)./2;

height_vct = tand(angles_of_base) * base_half_width;

z_vct_full = [];
for ii = 1 : length(angles)
    z_vct_half_tmp = tand(angles_of_base(ii)) * [0 : base_half_width];
    z_vct_full_tmp = [z_vct_half_tmp z_vct_half_tmp(end-1 : -1 :1)];
    z_vct_full = [z_vct_full z_vct_full_tmp(1 : end-1)]+1; 
end;

% plot(z_vct_full)
Z_to_test_angle_im = repmat(z_vct_full(:), 1, 176);
Z_to_test_angle_im=round(Z_to_test_angle_im);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_half_width = 8 ;
angles = [10, 20, 30, 40, 50, 60, 70, 80, 90];

[Z_to_test_angle_im, z_vct_full] = create_testing_zMaps(base_half_width, angles);

% figure; image(Z_to_test_angle_im); title('the depth data')
% colormap(gray(255))
% size(Z_to_test_angle_im);
% figure; plot(z_vct_full); title('the proflle of the depth, note the angles');

%%%%%%% add noise to the CG image %%%%%%%%%
%// Adjust intensities in image I to range from 0 to 1
im = Z_to_test_angle_im - min(Z_to_test_angle_im(:));
im = Z_to_test_angle_im / max(Z_to_test_angle_im(:));
% or
im=Z_to_test_angle_im/255;

%// Add noise to image
SNR=10;
v = var(im(:)) / 10^(SNR/10) ; 
%SNR = 10log10[var(image)/var(noise)]
%For a given image and SNR=5db, the variance of the noise would be:
%var(noise) = var(image)/10^(SNR/10) = var(image)/sqrt(10)
im_noisy = imnoise(im, 'gaussian', 0, v);
imwrite(im_noisy,'saw_tooth_SNR=10.bmp');
figure;  imshow(im_noisy);
%// Show images

subplot(1, 2, 1), imshow(im), title('Original image')
subplot(1, 2, 2), imshow(im_noisy), title('Noisy image, SNR=5db')