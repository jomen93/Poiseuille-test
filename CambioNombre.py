import numpy as np 

k = raw_input()
x = np.loadtxt("x")
y = np.loadtxt("y")
u = np.loadtxt("u")
v = np.loadtxt("v")
rho = np.loadtxt("rho")
np.savetxt("x"+str(k)+".dat", x)
np.savetxt("y"+str(k)+".dat", y)
np.savetxt("u"+str(k)+".dat", u)
np.savetxt("v"+str(k)+".dat", v)
np.savetxt("rho"+str(k)+".dat", rho)

