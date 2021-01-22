import numpy as np
import pyefd

img_1 = np.array(
    [
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            191,
            64,
            127,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            0,
            0,
            0,
            127,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            64,
            0,
            0,
            0,
            0,
            64,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            191,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            64,
            127,
            64,
            64,
            0,
            0,
            64,
            191,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            191,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            127,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            64,
            0,
            0,
            127,
            255,
            255,
            191,
            64,
            0,
            0,
            0,
            0,
            0,
            64,
            127,
            127,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            191,
            0,
            0,
            0,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            191,
            0,
            0,
            0,
            64,
            127,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            64,
            0,
            0,
            0,
            0,
            0,
            64,
            191,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            127,
            64,
            0,
            0,
            0,
            0,
            64,
            191,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            191,
            127,
            0,
            0,
            0,
            0,
            127,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            191,
            127,
            0,
            0,
            0,
            64,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            0,
            0,
            0,
            191,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            127,
            0,
            0,
            127,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            127,
            0,
            0,
            127,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            127,
            191,
            255,
            255,
            255,
            255,
            127,
            0,
            0,
            0,
            191,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            127,
            0,
            127,
            255,
            255,
            191,
            64,
            0,
            0,
            0,
            191,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            191,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            127,
            0,
            0,
            0,
            0,
            0,
            0,
            64,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            127,
            0,
            0,
            0,
            64,
            191,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
        [
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
            255,
        ],
    ]
)

import cv2

img_1 = np.uint8(img_1)

cv2.imwrite("test_raw.png", img_1)

edges = cv2.Canny(img_1, 100, 200)
contour_2 = []

for i in range(edges.shape[0]):
    for j in range(edges.shape[1]):
        if edges[i, j] == 255:
            contour_2.append([i, j])
contour_2 = np.array(contour_2)

cv2.imwrite("test1.png", img_1)

coeffs = pyefd.elliptic_fourier_descriptors(contour_2, order=10, normalize=False)

contour_2 = pyefd.reconstruct_contour(coeffs, locus=(0, 0), num_points=300)

for i in range(contour_1.shape[0]):
    tmp[int(round(contour_1[i][0]))][int(round(contour_1[i][1]))] = 255
print(tmp.shape)
cv2.imwrite("test2.png", tmp)