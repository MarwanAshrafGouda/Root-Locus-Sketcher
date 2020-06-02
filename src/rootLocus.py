import matplotlib.pyplot as plt
import numpy as np
import math
import sympy as sp

plt.plot([0, 0], [-80, 80], 'k-', lineWidth=0.6)
plt.plot([-100, 50], [0, 0], 'k-', lineWidth=0.6)

# plotting poles
poleX = [0, -25, -50, -50]
poleY = [0, 0, 10, -10]
plt.plot(poleX, poleY, 'bx', label='Pole')

# drawing asymptotic lines
asymptoteCount = len(poleX)
sumPoles = 0
angles = []
for i in range(0, asymptoteCount):
    sumPoles += complex(poleX[i], poleY[i])
    angle = (2 * i + 1) * math.pi / asymptoteCount
    angles.append(angle)
print("Asymptote angles:")
print(angles)
centroid = np.real(sumPoles / asymptoteCount)
plt.plot(centroid, 0, 'ro', label='Centroid')
print("Centroid:")
print(centroid)
for angle in angles:
    endX = 100 * math.cos(angle) + centroid
    endY = 100 * math.sin(angle)
    asymX = [centroid, endX]
    asymY = [0, endY]
    plt.plot(asymX, asymY, '--', lineWidth=1, label='Asymptote')

# finding break away points
s = sp.symbols('s')
exp = 1
for i in range(0, len(poleX)):
    exp *= (s - complex(poleX[i], poleY[i]))
exp = sp.expand(exp)
print("Characteristic Equation:")
print(exp)
expDiff = sp.Derivative(exp, s).doit()
print("Derivative of the Characteristic Equation")
print(expDiff)
roots = sp.solve(expDiff, s)
print("Roots of the Derivative of the Characteristic Equation")
print(roots)
realPoles = []
for i in range(0, len(poleY)):
    if poleY[i] == 0:
        realPoles.append(poleX[i])
realPoles = sorted(realPoles, reverse=True)
breakAwayPoints = []
for i in range(0, len(realPoles), 2):
    for root in roots:
        root = complex(root)
        if round(root.imag, 16) == 0:
            if realPoles[i] > root.real > realPoles[i + 1]:
                breakAwayPoints.append(root.real)
print("Break Away Points:")
print(breakAwayPoints)
plt.plot(breakAwayPoints, 0, 'm.', label='Break Away Point')

# finding imaginary axis intersection points
routhTable = []
for i in range(-1, len(poleX)):
    routhRow = []
    routhTable.append(routhRow)
for i in range(len(poleX), 0, -2):
    routhTable[0].append(exp.coeff(s, i))
for i in range(len(poleX) - 1, 0, -2):
    routhTable[1].append(exp.coeff(s, i))
k = sp.symbols('k')
routhTable[0].append(k)
routhTable[1].append(0)

print(routhTable)
plt.grid(True)
plt.legend(loc="lower center", ncol=3)
plt.show()
