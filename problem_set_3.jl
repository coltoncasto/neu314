############################################################################
# Problem Set 3
# NEU 314
# Author: Colton Casto
#
############################################################################

using PyPlot; using JLD

# Problem 1
## part a
image = imread("el-capitan.png")
imshow(image)

## part b writing a function
function image_colors(filename)
    # load and display image
    image = imread(filename)
    imshow(image)

    # extract each color coponent of image
    red   = image[:,:,1]
    green = image[:,:,2]
    blue  = image[:,:,3]

    return red, green, blue
end

## part e
red, green, blue = image_colors("el-capitan.png")
image2 = zeros(360, 640, 3)
image2[:,:,1] = green
image2[:,:,2] = blue
image2[:,:,3] = red

# plot image and image2 next to one another
figure(figsize = [10, 5])
title("RGB vs GBR Image")
subplot(1,2,1)
imshow(image)
subplot(1,2,2)
imshow(image2)

 # Problem 2
 ##
