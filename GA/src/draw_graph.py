import numpy as np
import matplotlib.pyplot as plt

scores = []
with open("sin.txt") as f:
    for line in f.readlines():
        line = line.strip()
        length = len(line)
        print(length, 2**length, int(line, 2), line)
        scores.append((int(line, 2) / (2**length)) * 128)
scores = np.array(scores)
print(scores)
x = np.arange(0, 128, 0.1)
y = np.sin(x)
plt.plot(x, y)
plt.scatter(scores, np.sin(scores), color='orange')
plt.show()
