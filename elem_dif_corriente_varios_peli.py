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
w=500#Hz
mu=1.25663706*10**(-6)  #Permeabilidad del vacío (H*m-1)
c=299792458 #Velocidad de la luz (m/s)


b=w/(c*1.)
londa=c/(w*1.)
P=1./w      #Periodo

#Parámetros respecto los otros elementos
des=np.pi/2  #Desfase de la intensidad
paso=4          #Paso de longitud de onda entre los elementos
dist=londa/(paso*1.)   #Distancia entre los elementos
num_elem=3      #Numero de elementos

#Discretizar el espacio
esplim=londa*13 #Límites del espacio. Cuando r mucho mayor que longitud onda
espp=esplim/150.  #Paso de Discretización del espacio. 300x300 cuadrados
cent=1  #Cuánto eliminar de la zona central de cada elemento
grafmax=0.2

dy=np.arange(-esplim,esplim,espp)
dz=np.arange(-esplim,esplim,espp)
cy,cz=np.meshgrid(dy,dz)

#Discretizar el tiempo de tal forma que en la película se vean al menos 5 periodos

t=np.arange(0,19*P/2.,P/7.)   #En cada periodod 10 valores diferentes del tiempo

#Calcular el radio y angulo a todos los elementos
r=np.zeros([num_elem,np.shape(cy)[0],np.shape(cy)[1]])
ang=np.zeros([num_elem,np.shape(cy)[0],np.shape(cy)[1]])
for i in range(0,num_elem,1):
    r[i,:]=np.sqrt((cy-dist*i)**2+cz**2)
    ang[i,:]=np.arctan2(cy-dist*i,cz)   #Se resta porque la distancia (en el eje y) al punto es menor si está a la izq


# arctan2 da el ángulo entre -pi y pi. Se pasa de 0 a 2pi
for index,j in np.ndenumerate(ang):
    if j<0:
        ang[index]=ang[index]+2*np.pi
        



##############################################################################
##                                                                          ##
##                           FUNCIONES                                      ##
##                                                                          ##
##############################################################################

##############################################################################
#Función que calcula el campo electrico

def E(t,r,ang,I,z,w,des):
    mu=1.25663706*10**(-6)  #Permeabilidad del vacío (H*m-1)
    c=299792458 #Velocidad de la luz (m/s)
    b=w/(c*1.)
# Solo el campo de radiacion (y por tanto solo componente en ang) y en el dominio del tiempo
    Elec=np.real((I*z*1j*w*mu*np.exp(-1j*b*r)*np.exp(1j*w*t)*np.sin(ang)*np.exp(1j*des))/(4*np.pi*r))
    return Elec    


def E_elem(t,r,ang,I,z,w,num_elem,des):
    #Calcular un primer elemento para no tenrer q hacer una matriz nueva en cada tiempo
    E_tot=E(t,r[0,:,:],ang[0,:,:],I,z,w,des*0)
    for i in range(1,num_elem,1):
        E_elem=E(t,r[i,:,:],ang[i,:,:],I,z,w,des*i)
        E_tot=E_tot+E_elem
        
        #Quitar el centro de cada elemento situado a i*londa/paso. Solo se mira la primera fila de cy porque está ordenado por columnas
        #No se puede meter en el anterior loop porque si no suma encima del cero
        cent_pos=np.argwhere(((-espp*cent)<dy) & (dy<(+espp*cent)))
        E_tot[np.shape(dz)[0]/2-cent-1:np.shape(dz)[0]/2+cent+1,cent_pos]=0     #Se quitan las posiciones calculadas en y y un numero fijo en z

    for i in range(1,num_elem,1):    
        cent_pos=np.argwhere(((dist*i-espp*cent)<=dy) & (dy<=(dist*i+espp*cent)))
        E_tot[np.shape(dz)[0]/2-cent-1:np.shape(dz)[0]/2+cent+1,cent_pos]=0     #Se quitan las posiciones calculadas en y y un numero fijo en z
    
    return np.abs(E_tot) #Tomamos el valor absoluto porque el signo solo indica la dirección del vector theta 
    
#############################################################################
# Función que grafica el campo eléctrico para cada tiempo en un mapa de color
# Es necesaria para plotearlo variable en el tiempo
def GrafEColor(t,r,ang,I,z,w,Emax,num_elem,des):
    
    Etot=E_elem(t,r,ang,I,z,w,num_elem,des)
    Etot=Etot/(Emax*1.)  #Se normaliza al un valor máximo inicial
    Etot[Etot>1]=1

    graf=ax.pcolormesh(dy/(londa*1.),dz/(londa*1.),Etot,vmin=0,vmax=grafmax,cmap='inferno')  #vmax para ver bien los colores
    ax.set_title('E-field %d elem., distancia=$\lambda$/%d , $\phi$=%.3f [rad]\n $\omega$=%.f [Hz]; $\lambda$=%.2e [m]; t=%.2e [s]' %(num_elem,paso,des,w,londa,t))
    return graf
    
#############################################################################
# Función que grafica el campo eléctrico para cada tiempo en 3D
# Es necesaria para plotearlo variable en el tiempo
def GrafE3D(t,r,ang,I,z,w,Emax,num_elem,des):
    
    Etot=E_elem(t,r,ang,I,z,w,num_elem,des) 
    Etot=Etot/(Emax*1.)  #Se normaliza al un valor máximo inicial    
    Etot[Etot>1]=1

    ax.clear()   
    ax.set_xlabel("Eje y [y/$\lambda$]")
    ax.set_ylabel("Eje z [z/$\lambda$]")
    ax.set_zlabel("E/$E_{0max}$")
    ax.set_zlim(0,1)

    ax.view_init(elev=50)
    graf=ax.plot_surface(cy/(londa*1.),cz/(londa*1.),Etot,vmax=1)
    ax.set_title('E-field %d elem., distancia=$\lambda$/%d, $\phi$=%.3f [rad]\n $\omega$=%.f [Hz]; $\lambda$=%.2e [m]; t=%.2e [s]' %(num_elem,paso,des,w,londa,t))

    return graf

##############################################################################
##                                                                          ##
##                  PROGRAMA PRINCIPAL                                      ##
##                                                                          ##
##############################################################################

#Se calcula un primer término para hallar el máx que la gráfica debe representar para que no cambie
Emax=E_elem(t[0],r,ang,I,z,w,num_elem,des)
Emax=np.max(Emax)

'''
#######################################################
#Un solo gráfico
#######################################################

#Se crea un entorno gráfico inicial

fig,ax=plt.subplots(figsize=(50,100))

ax.set_xlabel("Eje y [y/$\lambda$]")
ax.set_ylabel("Eje z [z/$\lambda$]")

GrafEColor(0.05,r,ang,I,z,w,Emax,num_elem,des)
plt.savefig('../Elemento_dif_corriente/Varios_elementos/elem%d_londa%d_des%d_grafica.mp4' %(num_elem,paso,des*100))
'''

#######################################################
#Gráfico con campo de color
#######################################################

#Se crea un entorno gráfico inicial

fig,ax=plt.subplots(figsize=(50,100))

ax.set_xlabel("Eje y [y/$\lambda$]")
ax.set_ylabel("Eje z [z/$\lambda$]")

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1,codec="libx264")

#Se dibuja para cada tiempo la gráfica
peli=animation.FuncAnimation(fig,GrafEColor,fargs=(r,ang,I,z,w,Emax,num_elem,des),frames=t,interval=10,repeat=False)

FFwriter = animation.FFMpegWriter(fps=15,extra_args=['-vcodec', 'libx264'])
peli.save('../Elemento_dif_corriente/Varios_elementos/elem%d_londa%d_des%d_color.mp4' %(num_elem,paso,des*100), writer = FFwriter)

print "Película terminada\n\n"

'''
#######################################################
#Gráfico en 3D
#######################################################

#Se crea un entorno gráfico inicial
fig = plt.figure(figsize=(50,100))
ax = Axes3D(fig)

cent=cent     #En este caso se ve mejor si se elimina más parte central del elemento

#Formato de la película
Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1,codec="libx264")

#Se dibuja para cada tiempo la gráfica
peli=animation.FuncAnimation(fig,GrafE3D,fargs=(r,ang,I,z,w,Emax,num_elem,des),frames=t,interval=10,repeat=False)

FFwriter = animation.FFMpegWriter(fps=15,extra_args=['-vcodec', 'libx264'])
peli.save('../Elemento_dif_corriente/Varios_elementos/elem%d_londa%d_des%d_3d.mp4' %(num_elem,paso,des*100), writer = FFwriter)

print "Película terminada\n\n"
'''
