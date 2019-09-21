# u1, u2, sigma = -3, 3, 1
u1, u2, sigma = -3, 3, 0.1

w1 = (u1-u2)/sigma**2
w2 = (u2**2 - u1**2)/2*sigma**2
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-.005, 0.005, 1000)

a = w1*x + w2

f = 1/(1 + np.exp(-a))
fig = plt.figure(figsize=(6, 5))
plt.title(f'm1 = {u1}, m2 = {u2} and sigma = {sigma}', fontsize = 15)
plt.plot(x, a, label=f'a(x) = {round(w1)}x + {w2}')
plt.plot(x, f, label='P(z = 1 | x, theta)')
plt.ylabel('y', fontsize = 15)
plt.xlabel('x', fontsize = 15)
plt.legend(fontsize = 12)
plt.show()
fig.savefig(f'sigma_01.png', bbox_inches='tight')