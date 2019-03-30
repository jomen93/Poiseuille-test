import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm 


def uy(x,P):
    d =1.0; eta = 17/60.
    return (1/(2*eta))*P*((x-d/2.)**2 - d**2/4.)

# grafica para graficar los datos obtenidos desde c++

def uy(x,P):
    d =1.0; eta = 17/60.
    return (1/(2*eta))*P*((x-d/2.)**2 - d**2/4.)

def PoiseuilleGrafica(udat,vdat,xdat,ydat,rhodat,pasos):
    
    u = np.transpose(np.loadtxt(udat, unpack = True))
    v = np.transpose(np.loadtxt(vdat, unpack = True))
    x = np.transpose(np.loadtxt(xdat, unpack = True))
    y = np.transpose(np.loadtxt(ydat, unpack = True))
    rho = np.transpose(np.loadtxt(rhodat, unpack = True))
    X = np.linspace(0,1,50)


    f, axarr = plt.subplots(1,2, figsize=(20,4))
    st = f.suptitle("Flujo de Pouseuille $\\tau = 0.9$", fontsize="x-large")
    st.set_y(0.95)

    M = np.hypot(u, v)
    axarr[0].streamplot(x,y,u,v, color="k",linewidth=0.8,density=0.4, arrowstyle='->', arrowsize=1.5)
    axarr[0].quiver(x, y, u,v, M , cmap=plt.cm.jet, width=0.022,scale=1/0.1)
    axarr[0].set_title(str(pasos)+" pasos")
    #axarr[0].set_xticklabels([])
    
    axarr[1].plot(x,v[128,:],"b", label = "Simulacion")
    #axarr[1].plot(x,y[128,:],"b", label = "Simulacion")
    if (k == 100000):
    	axarr[1].plot(X,uy(X,0.65),"r+", label = "Teorica")
    axarr[1].legend()
    axarr[1].grid(True)
    axarr[1].set_title('Perfil de Velocidad')

    im = axarr[0].quiver(x,y, u,v, M , cmap=plt.cm.jet, width=0.022,scale=1/0.1)
    f.colorbar(im)#, ax=axarr, shrink = 1.0)
    plt.savefig(str(k)+"_pasos")

k = raw_input()
PoiseuilleGrafica("Datos/u"+str(k)+".dat","Datos/v"+str(k)+".dat","Datos/x"+str(k)+".dat","Datos/y"+str(k)+".dat","Datos/rho"+str(k)+".dat",k)
