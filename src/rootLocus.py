import matplotlib.pyplot as plt
import numpy as np
import math
import sympy as sp

plt.plot([0, 0], [-60, 60], 'k-', lineWidth=0.6)
plt.plot([-90, 20], [0, 0], 'k-', lineWidth=0.6)

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
    endX = 72 * math.cos(angle) + centroid
    endY = 72 * math.sin(angle)
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
for i in range(1, len(routhTable) - 1):
    for j in range(0, len(routhTable[i])):
        if j == len(routhTable[i]) - 1:
            routhTable[i + 1].append(0)
        else:
            routhTable[i + 1].append(
                (routhTable[i][j] * routhTable[i - 1][j + 1] - routhTable[i][j + 1] * routhTable[i - 1][j]) /
                routhTable[i][
                    j])
print("Routh Table:")
print(routhTable)
kVal = sp.solve(routhTable[len(routhTable) - 2])
print("Critical Value of K:")
criticalK = kVal[k]
print(criticalK)
auxEqu = routhTable[len(routhTable) - 3][0] * s ** 2 + criticalK
print("Auxiliary Equation:")
print(auxEqu)
imaginaryAxisIntercepts = sp.solve(auxEqu, s)
print("Imaginary Axis Intercepts:")
print(imaginaryAxisIntercepts)
for intercept in imaginaryAxisIntercepts:
    intercept = complex(intercept)
    plt.plot(0, intercept.imag, 'c.', label='Imaginary Axis Intersection Point')

# calculating angle of departure
departureAngles = []
for i in range(0, len(poleY)):
    sigmaPhi = 0
    for j in range(0, len(poleY)):
        if i != j:
            phi = complex(poleX[i], poleY[i]) - complex(poleX[j], poleY[j])
            if phi.real == 0:
                if phi.imag > 0:
                    sigmaPhi += math.pi / 2
                if phi.imag < 0:
                    sigmaPhi += 3 * math.pi / 2
            elif phi.imag == 0:
                if phi.real > 0:
                    sigmaPhi += 0
                if phi.real < 0:
                    sigmaPhi += math.pi
            else:
                sigmaPhi += math.atan(phi.imag / phi.real)
    departureAngle = math.pi - sigmaPhi
    if departureAngle < 0:
        departureAngle += 2 * math.pi
    departureAngles.append(departureAngle)
print("Departure Angles:")
print(departureAngles)
radius = 3
for i in range(0, len(departureAngles)):
    angleX = []
    angleY = []
    theta = 0
    while theta < departureAngles[i]:
        angleX.append(radius * math.cos(theta) + poleX[i])
        angleY.append(radius * math.sin(theta) + poleY[i])
        theta += departureAngles[i] / 100
    plt.plot(angleX, angleY, 'y', label='Angle of Departure', lineWidth=1)

# drawing root locus
rootsX = []
rootsY = []
print("This will take a while; please wait...")
flippingPoint = -exp.subs(s, breakAwayPoints[0])
for k in range(0, 10000000, 2000):
    roots = sp.solve(exp + k, s, simplify=False, rational=False)
    rootX = []
    rootY = []
    for root in roots:
        rootX.append(complex(root).real)
        rootY.append(complex(root).imag)
    if k > flippingPoint:
        rootX.append(rootX[0])
        rootX.append(rootX[1])
        del rootX[0]
        del rootX[0]
        rootY.append(rootY[0])
        rootY.append(rootY[1])
        del rootY[0]
        del rootY[0]
    rootsX.append(rootX)
    rootsY.append(rootY)
for i in range(0, len(rootsX[0])):
    colX = [row[i] for row in rootsX]
    colY = [row[i] for row in rootsY]
    plt.plot(colX, colY, label='Root Locus')

plt.grid(True)
plt.legend(loc="lower center", ncol=3)
plt.show()
