############################################################################
# Problem Set 3
# NEU 314
# Author: Colton Casto
#
############################################################################

using PyPlot; using JLD

# Problem 1 part 1
## part a
image = imread("el-capitan.png")
imshow(image)

## part b writing a function
"""
image_colors() -- extracts color components of an image.
This function loads an image given a file name and extracts the RGB components
    of the image. The components are then returned individually.

Args:
   filename (str): name of file where the image is stored

Returns:
   red (array): red component of original image with dimesnions row x col
   green (array): green component of original image with dimesnions row x col
   blue (array): blue component of original image with dimesnions row x col
"""
function image_colors(filename)
    # load and display image
    println("Reset with git reset --mixed")
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

 ## Problem 1 part 2
 """
circular_rotation() -- circularly moves one channel (red channel) up p number
                       of pixeles.
This function takes an image array and and integer p and circularly rotates
    one chanel up some number of pixeles, such that the top p rows now become
    the bottom p rows.

Args:
   image (array): Image array with dimensions row x col x color(and alpha).
                  Color order is RGB.
   p (int): number of rows to rotate circularly

Returns:
   out (array): Image array with same dimensions as input image but with
                the red channel circularly rotated p rows.
"""
 function circular_rotation(image, p)
     nrows = size(image)[1]
     circ_image = copy(image)
     circ_image[1:(nrows-p),:,1] = image[(p+1):nrows,:,1]
     circ_image[(nrows-p+1):nrows,:,1] = image[1:p,:,1]
     return circ_image
 end

 circular_rotation(image, 100)
