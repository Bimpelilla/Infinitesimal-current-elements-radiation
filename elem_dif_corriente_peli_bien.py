# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 16:36:16 2017

@author: M. Pedrosa Bustos
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


##############################################################################
##                                                                          ##
##                     DEFINICIÓN DE PARÁMETROS                             ##
##                                                                          ##
##############################################################################

I=100
z=1
w=50#Hz
mu=1.25663706*10**(-6)  #Permeabilidad del vacío (H*m-1)
c=299792458 #Velocidad de la luz (m/s)

b=w/(c*1.)
londa=c/(w*1.)

#Discretizar el espacio
esplim=londa*17 #Límites del espacio. Cuando r mucho mayor que longitud onda
espp=esplim/150.  #Paso de Discretización del espacio. 300x300 cuadrados
cent=0  #Eliminar zona central

dy=np.arange(-esplim,esplim,espp)
dz=np.arange(-esplim,esplim,espp)
cy,cz=np.meshgrid(dy,dz)
r=np.sqrt(cy**2+cz**2)
ang=np.arctan2(cy,cz)

# arctan2 da el ángulo entre -pi y pi. Se pasa de 0 a 2pi
for index,i in np.ndenumerate(ang):
    if i<0:
        ang[index]=ang[index]+2*np.pi
        
#Discretizar el tiempo de tal forma que en la película se vean al menos 5 periodos
P=1./w      #Periodo
t=np.arange(0,10*P,P/7.)   #En cada periodod 10 valores diferentes del tiempo



##############################################################################
##                                                                          ##
##                           FUNCIONES                                      ##
##                                                                          ##
##############################################################################

##############################################################################
#Función que calcula el campo electrico

def E(t,r,ang,I,z,w):
    mu=1.25663706*10**(-6)  #Permeabilidad del vacío (H*m-1)
    c=299792458 #Velocidad de la luz (m/s)
    b=w/(c*1.)
# Solo el campo de radiacion (y por tanto solo componente en ang) y en el dominio del tiempo
    Elec=np.real((I*z*1j*w*mu*np.exp(-1j*b*r)*np.exp(1j*w*t)*np.sin(ang))/(4*np.pi*r))
    Elec=np.abs(Elec)   #Tomamos parte real y el valor absoluto porque el signo solo indica la dirección del vector theta 
    
    #Se quita la zona central donde esta el elemento de corriente para que no explote
    cent_pos=np.argwhere(((0.-espp*cent-1)<=dy) & (dy<=(0.+espp*cent)))
    Elec[np.shape(dz)[0]/2-cent-1:np.shape(dz)[0]/2+cent,cent_pos]=0     #Se quitan las posiciones calculadas en y y un numero fijo en z

    return Elec    


#############################################################################
# Función que grafica el campo eléctrico para cada tiempo en un mapa de color
# Es necesaria para plotearlo variable en el tiempo
def GrafEColor(t,r,ang,I,z,w,Emax):
    
    Etot=E(t,r,ang,I,z,w) 
    Etot=Etot/(Emax*1.)  #Se normaliza al un valor máximo inicial
    Etot[Etot>1.]=1

    graf=ax.pcolormesh(dy/(londa*1.),dz/(londa*1.),Etot,vmin=0,vmax=1,cmap='inferno')  #vmax para ver bien los colores
    ax.set_title('E field\n $\omega$=%.f [Hz]; $\lambda$=%.2e [m]; t=%.2e [s]' %(w,londa,t))
    return graf
    
#############################################################################
# Función que grafica el campo eléctrico para cada tiempo en 3D
# Es necesaria para plotearlo variable en el tiempo
def GrafE3D(t,r,ang,I,z,w,Emax):
    
    Etot=E(t,r,ang,I,z,w) 
    Etot=Etot/(Emax*1.)  #Se normaliza al un valor máximo inicial
    Etot[Etot>1.]=1    
    ax.clear()   
    ax.set_xlabel("Eje y [y/$\lambda$]")
    ax.set_ylabel("Eje z [z/$\lambda$]")
    ax.set_zlabel("E/$E_{0max}$")
    ax.set_zlim(0,1)

    ax.view_init(elev=50)
    graf=ax.plot_surface(cy/(londa*1.),cz/(londa*1.),Etot)
    ax.set_title('E field\n $\omega$=%.f [Hz]; $\lambda$=%.2e [m]; t=%.2e [s]' %(w,londa,t))

    return graf

##############################################################################
##                                                                          ##
##                  PROGRAMA PRINCIPAL                                      ##
##                                                                          ##
##############################################################################

#Se calcula un primer término para hallar el máx que la gráfica debe representar para que no cambie
#Se coge el maximo de la mitad del periodo
Emax=E(t[0],r,ang,I,z,w)

Emax=np.max(Emax)


#######################################################
#Gráfico con campo de color
#######################################################

#Se crea un entorno gráfico inicial

fig,ax=plt.subplots(figsize=(50,100))

ax.set_xlabel("Eje y [y/$\lambda$]")
ax.set_ylabel("Eje z [z/$\lambda$]")


# Formato para la película
Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1,codec="libx264")

#Se dibuja para cada tiempo la gráfica
peli=animation.FuncAnimation(fig,GrafEColor,fargs=(r,ang,I,z,w,Emax),frames=t,interval=10,repeat=False)

FFwriter = animation.FFMpegWriter(fps=15,extra_args=['-vcodec', 'libx264'])
peli.save('../Elemento_dif_corriente/1_elemento/1elem_color.mp4', writer = FFwriter)

print "Película terminada\n\n"

'''
#######################################################
#Gráfico en 3D
#######################################################

#Se crea un entorno gráfico inicial
fig = plt.figure(figsize=(50,100))
ax = Axes3D(fig)
cent=cent+7 #Se quita más del centro para que se vea mejor

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1,codec="libx264")

#Se dibuja para cada tiempo la gráfica
peli=animation.FuncAnimation(fig,GrafE3D,fargs=(r,ang,I,z,w,Emax),frames=t,interval=10,repeat=False)

FFwriter = animation.FFMpegWriter(fps=15,extra_args=['-vcodec', 'libx264'])
peli.save('../Elemento_dif_corriente/1_elemento/1elem_3d.mp4', writer = FFwriter)

print "Película terminada\n\n"
'''
