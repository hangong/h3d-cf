# h3d-cf

The code is provided for the following publications: 

[1] Gong, H., Finlayson, G.D., Fisher, R.B. and Fang, F., 2017. 3D color homography model for photo-realistic color transfer re-coding. The Visual Computer, pp.1-11.

[2] Gong, H., Finlayson, G.D., Fisher, R.B.: Recoding color transfer as a color homography. In: British Machine Vision Conference. BMVA (2016)

Please cite the above papers if you find them useful for your research.

Copyright 2018 Han Gong (<gong@fedoraproject.org>), University of {East Anglia, Edinburgh}.

# Introduction
Color transfer is an image editing process that naturally transfers the color theme of a source image to a target image. In this paper, we propose a 3D color homography model which approximates photo-realistic color transfer algorithm as a combination of a 3D perspective transform and a mean intensity mapping. A key advantage of our approach is that the re-coded color transfer algorithm is simple and accurate. Our evaluation demonstrates that our 3D color homography model delivers leading color transfer re-coding performance. In addition, we also show that our 3D color homography model can be applied to color transfer artifact fixing, complex color transfer acceleration, and color-robust image stitching.

![pipeline](https://media.springernature.com/original/springer-static/image/art%3A10.1007%2Fs00371-017-1462-x/MediaObjects/371_2017_1462_Fig1_HTML.gif)

# Usage

* `cf_3D_H.m`: the code of the TVC paper [1].
* `cf_2D_H.m`: the code of the BMVC paper [2].

# Testing
* `test_main.m`: this file contains the evaluation code for re-producing the results in the paper.

`cf_MK.m` and `cf_Poly.m` are the implementations for two other color transfer approximation methods for comparison. The approximation outputs will be saved in the folder `out_pair`.

# Dataset
A set of 200 transfer cases for 4 different color transfer methods (800 test cases in total) are contained in the folder `in_pair`. Note that the code for the *original* color transfer methods is not provided here. Only the cached results for color transfer approximation is provided in the dataset.

Rules of Filenames:
* `*_s.jpg`: original source image.
* `*_s.jpg`: original target image.
* `*_[method_name].jpg`: color transfer result produced by [method_name].

# Applications
The proposed method may be used for (see the original paper for higher-quality images):
* color transfer artifact fixing
![pipeline](https://media.springernature.com/lw785/springer-static/image/art%3A10.1007%2Fs00371-017-1462-x/MediaObjects/371_2017_1462_Fig6_HTML.gif)
* complex color transfer acceleration
![pipeline](https://media.springernature.com/lw785/springer-static/image/art%3A10.1007%2Fs00371-017-1462-x/MediaObjects/371_2017_1462_Fig5_HTML.gif)
* color-robust image stitching.
![pipeline](https://media.springernature.com/lw785/springer-static/image/art%3A10.1007%2Fs00371-017-1462-x/MediaObjects/371_2017_1462_Fig7_HTML.gif)
