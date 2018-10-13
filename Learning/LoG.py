import cv2 as cv
import numpy as np
from matplotlib import pyplot as plt

image = cv.imread('FM15-M13GFP-2.png', 0)

image_edge_detect = cv.GaussianBlur(image,(5,5), 0.8)
laplacian = cv.Laplacian(image_edge_detect , cv.CV_64F,  ksize=31)# ksize should be an odd
La = np.array(laplacian)
normalize = cv.normalize(laplacian,laplacian,0,255,cv.NORM_MINMAX)
no = np.array(normalize)
mor = cv.morphologyEx(normalize, cv.MORPH_GRADIENT, np.ones((5,5),np.uint8), iterations= 2)
#mor = cv.morphologyEx(normalize, cv.MORPH_CLOSE, np.ones((5,5),np.uint8), iterations= 2)

#mor = cv.morphologyEx(mor, cv.MORPH_OPEN, np.ones((5,5),np.uint8), iterations= 2)
#the = cv.adaptiveThreshold(mor, 255, cv.ADAPTIVE_THRESH_GAUSSIAN_C,\
#                           cv.THRESH_BINARY, 11, 2)
#edged, contours, hierarchy = cv.findContours(mor, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)

#filteredcnt = [ i for i in contours if cv.contourArea(i)>50]
#edge = cv.drawContours(image.copy(),filteredcnt,-1,(0,255,0),3)
#edges = cv.drawContours(image.copy(), filteredcnt, -1, (0,255,255),3)
#print(np.array(filteredcnt).shape)U


plt.figure(figsize=(16,16))
plt.subplot(221)
plt.imshow(image, cmap = 'jet')
plt.subplot(222)
plt.imshow(laplacian, cmap='gray')
plt.subplot(223)
plt.imshow(normalize, cmap='gray')
plt.subplot(224)
plt.imshow(mor, cmap='gray')
plt.show()