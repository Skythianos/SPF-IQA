# No-Reference Image Quality Assessment Based on the Fusion of Statistical and Perceptual Features
The goal of no-reference image quality assessment (NR-IQA) is to predict the quality of an image as perceived by human observers without using any pristine, reference images. In this study, an NR-IQA algorithm is proposed which is driven by a novel feature vector containing statistical and perceptual features. Different from other methods, normalized local fractal dimension distribution and normalized first digit distributions in the wavelet and spatial domains are incorporated into the statistical features. Moreover, powerful perceptual features, such as colorfulness, dark channel feature, entropy, and mean of phase congruency image, are also incorporated to the proposed model. Experimental results on five large publicly available databases (KADID-10k, ESPL-LIVE HDR, CSIQ, TID2013, and TID2008) show that the proposed method is able to outperform other state-of-the-art methods.

If you use this MATLAB code please cite the following paper:<br/>
@article{varga2020no,<br/>
  title={No-Reference Image Quality Assessment Based on the Fusion of Statistical and Perceptual Features},<br/>
  author={Varga, Domonkos},<br/>
  journal={Journal of Imaging},<br/>
  volume={6},<br/>
  number={8},<br/>
  pages={75},<br/>
  year={2020},<br/>
  publisher={Multidisciplinary Digital Publishing Institute}<br/>
}

This code was written and tested in MATLAB R2019a. Required toolboxes:<br/>
- System Identification Toolbox<br/>
- Image Processing Toolbox<br/>
- Statistics and Machine Learning Toolbox<br/>
- Wavelet Toolbox<br/>
