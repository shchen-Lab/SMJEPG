import cv2
import numpy as np

with open('bird.txt') as fb:
  lines = fb.readlines()
  x = []
  for line in lines:
    x.append(int(line))

x = np.array(x).astype(np.uint8).reshape(1200, 900, 3)

cv2.imwrite('bird_test.jpg', x)
