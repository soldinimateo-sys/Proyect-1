import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Devuelve la fila (anterior) i donde se encuentran las propiedades termodinámicas correspondientes al nivel de energía (un) 
"""
def iteracion(un, v, i, data):
    x1 = (v-data[i,2])/(data[i,3]-data[i,2])
    u1 = data[i,4]+x1*(data[i,5]-data[i,4])
    """
    if un > u1:
        return iteracion(un,v,i-1,data)
    """
    x2 = (v-data[i+1,2])/(data[i+1,3]-data[i+1,2]) 
    u2 = data[i+1,4]+x2*(data[i+1,5]-data[i+1,4])
    while u2 > un:
        i+=1
        x2 = (v-data[i+1,2])/(data[i+1,3]-data[i+1,2]) 
        u2 = data[i+1,4]+x2*(data[i+1,5]-data[i+1,4])
    return i

"""
Dada una fila i, interpola la Temperatura o Presión correspondiente al nivel de enrgía (un)
"""
def interpol_PoT(j,i,un,v,data):
    "j = 0 para presión, j = 1 para Temperatura"
    x1 = (v-data[i,2])/(data[i,3]-data[i,2])
    u1 = data[i,4]+x1*(data[i,5]-data[i,4])
    x2 = (v-data[i+1,2])/(data[i+1,3]-data[i+1,2])
    u2 = data[i+1,4]+x2*(data[i+1,5]-data[i+1,4])
    x = [u2,u1]
    y = [data[i+1,j],data[i,j]]
    return np.interp(un,x,y)

"""
Función que determina la variación de Energía según la Temperatura a la que se encuentra el sistema
"""
def Qn(Tn, T0, data, dt, h_lat, h_sup = 0.006, r = 0.275, h = 0.9,):
    A_at = 2*np.pi*h*r
    A_sup = np.pi*(r**2)
    return -1*(h_lat*A_at+h_sup*A_sup)*(Tn-T0)*dt


# din definir: h_lat
def sim_Cooler(data, h_lat):
    dt = 5
    
    P = []
    T = []
    E = []

    # Ambiente:
    T0 = 20+273.15
    P0 = 1

    # Sistema de estudio:
    m = 10 # Masa de agua (m) en kilo gramos
    r = 0.55/2 # Radio Cilindro (r) en metros
    h = 0.9 # Altura Cilindro (h) en metros
    V = (r**2)*np.pi*h # Volumen Cilindro (V) en metros^3
    v = V/m

    # Estado Inicial
    i = 0
    tiempo = 0
    x = (v-data[i,2]) / (data[i,3]-data[i,2])
    u1 = data[i,4] + x * (data[i,5]-data[i,4])
    E.append(m * u1)
    T.append(data[i,1])
    P.append(data[i,0])

    # Pasan primeros 5 segundos ("loop anterior" al primero)
    tiempo += dt
    En = E[-1] + Qn(T[-1],T0,data, dt, h_lat)
    i = iteracion(En[-1]/m, v, i, data)
    Pn = interpol_PoT(0,i, En/m, v,data)

    while np.abs(Pn-P0) < 0.7:
        # Registro de (P,T,E) válidas
        P.append(Pn)
        T.append(interpol_PoT(1,i, En/m, v,data))
        E.append(En)
        # Pasan 5 segundos
        tiempo += dt
        En = E[-1] + Qn(T[-1],T0,data, dt, h_lat) 
        i = iteracion(En/m, v, i, data)
        Pn = interpol_PoT(0,i,En/m,v,data)

    return tiempo


