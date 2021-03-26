import numpy as np


x1 = np.array([[1, 2, 1]])
x2 = np.array([[1, 1, 1]])

# print(np.dot(x1, x2) / (np.sqrt(np.sum(x1 ** 2)) * np.sqrt(np.sum(x2 ** 2))))
print(np.dot(x1[0], x2[0]) / (np.sqrt(np.sum(x1 ** 2)) * np.sqrt(np.sum(x2 ** 2))))

theta_i = np.arccos(np.dot(x1[0], x2[0]) / (np.sqrt(np.sum(x1[0] ** 2)) * np.sqrt(np.sum(x2[0] ** 2))))

print(theta_i)