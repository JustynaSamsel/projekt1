# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:27:21 2022

@author: Justyna
"""



import numpy as np
from kodydoprojektu import *


# utworzenie instancji klasy Transformacje
elipsoida_grs80 = Transformacje(model = "grs80")

# nazwa pliku do odczytu
plik = "wsp_inp.txt"

# odczyt danych z pliku
tablica = np.genfromtxt(plik, delimiter=",", skip_header = 4)  # pomijam 4 pierwsze linie

w, k = np.shape(tablica)   # rozmiar tablicy (12,3)

# puste talice do obliczen
fi_lam_h = np.zeros((w, k))
XYZ = np.zeros((w,k))
XY_00 = np.zeros((w, 2))
XY_92 = np.zeros((w, 2))
n_e_u = np.zeros((w,k))
azymut_elewacja = np.zeros((11, 2))
odleglosc_2D_3D = np.zeros((11, 2))

# wyniki: (fi,lam,h), [x00, y00], [x92,y92] ---> 7 kolumn
tablica_wynikow = np.zeros((w, 7))


for ix in range(w):
    
    fi_lam_h[ix] = elipsoida_grs80.hirvonen(tablica[ix,0], tablica[ix,1], tablica[ix,2])
    XYZ[ix] = elipsoida_grs80.geodezyjne2XYZ(fi_lam_h[ix,0], fi_lam_h[ix,1], fi_lam_h[ix,2])
    XY_00[ix] = elipsoida_grs80.u2000(fi_lam_h[ix,0], fi_lam_h[ix,1])
    XY_92[ix] = elipsoida_grs80.u1992(fi_lam_h[ix,0], fi_lam_h[ix,1])
    n_e_u[ix] = elipsoida_grs80.neu(tablica[ix,0], tablica[ix,1], tablica[ix,2], tablica[ix,0]+1, tablica[ix,1]+1, tablica[ix,2]+1)
    tablica_wynikow[ix] = np.hstack([fi_lam_h[ix], XY_00[ix], XY_92[ix]])
# hstack --> automatyczne dopasowowywanie kolumn
 
for ix in range(w-1):
    odleglosc_2D_3D[ix] = elipsoida_grs80.odleglosc_2D_3D(tablica[ix,0], tablica[ix,1], tablica[ix,2], tablica[ix+1,0], tablica[ix+1,1], tablica[ix+1,2])
    azymut_elewacja[ix] = elipsoida_grs80.azymut_elewacja(tablica[ix,0], tablica[ix,1], tablica[ix,2], tablica[ix+1,0], tablica[ix+1,1], tablica[ix+1,2])  



### Panel dla uzytkownika

print("Wybierz jedną z poniższych opcji:")
while True:
    print("\nOpcje:")
    print("1 | Transformacja współrzędnych geocentrycznych na współrzędne geodezyjne")
    print("2 | Transformacja współrzędnych geodezyjnych na współrzędne geocentryczne")
    print("3 | Obliczenie współrzędnych topocentrycznych E,N,U")
    print("4 | Transformacja współrzędnych geodezyjnych do układu płaskiego 2000")
    print("5 | Transformacja współrzędnych geodezyjnych do układu płaskiego 1992")
    print("6 | Wyznaczenie kąta azymutu i kąta elewacji")
    print("7 | Wyznaczenie odległości 2D oraz 3D")

    opcja = input()

    if  opcja == "1":
        print('Szerokość geodezyjna | Długość geodezyjna | Wysokość na elipsoidzie')
        for ix in range(w):
            fi_lam_h[ix] = elipsoida_grs80.hirvonen(tablica[ix,0], tablica[ix,1], tablica[ix,2])
            print(f'{fi_lam_h[ix][0]} | {fi_lam_h[ix][1]} | {fi_lam_h[ix][2]}')
            
            
    elif  opcja == "2":
        print('X | Y | Z')
        for ix in range(w):
            XYZ[ix] = elipsoida_grs80.geodezyjne2XYZ(fi_lam_h[ix,0], fi_lam_h[ix,1], fi_lam_h[ix,2])
            print(f'{ XYZ[ix][0] } |{ XYZ[ix][1] }|{ XYZ[ix][2] }')
    
    
    elif  opcja == "3":
        print('N | E | U')
        for ix in range(w):
            n_e_u[ix] = elipsoida_grs80.neu(tablica[ix,0], tablica[ix,1], tablica[ix,2], tablica[ix,0]+1, tablica[ix,1]+1, tablica[ix,2]+1)
            print(f'{ n_e_u[ix][0] }|{ n_e_u[ix][1] }|{ n_e_u[ix][2] }')
            
            
            
    elif  opcja == "4":
        print('X_układ_2000|Y_układ_2000')
        for ix in range(w):
            XY_00[ix] = elipsoida_grs80.u2000(fi_lam_h[ix,0], fi_lam_h[ix,1])
            print(f'{ XY_00[ix][0] }|{ XY_00[ix][1] }')
    
    elif  opcja == "5":
        print('X_układ_1992|Y_układ_1992')
        for ix in range(w):
            XY_92[ix] = elipsoida_grs80.u1992(fi_lam_h[ix,0], fi_lam_h[ix,1])
            print(f'{ XY_92[ix][0] }|{ XY_92[ix][1] }')
            
            
    elif  opcja == "6":
        print('Azymut | Kąt elewacji')
        for ix in range(w-1):
            azymut_elewacja[ix] = elipsoida_grs80.azymut_elewacja(tablica[ix,0], tablica[ix,1], tablica[ix,2], tablica[ix+1,0], tablica[ix+1,1], tablica[ix+1,2])
            print(f'{ azymut_elewacja[ix][0] }|{ azymut_elewacja[ix][1] }')
    
    
    elif  opcja == "7":
        print('Odległość 2D | Odległość 3D')
        for ix in range(w-1):
            odleglosc_2D_3D[ix] = elipsoida_grs80.odleglosc_2D_3D(tablica[ix,0], tablica[ix,1], tablica[ix,2], tablica[ix+1,0], tablica[ix+1,1], tablica[ix+1,2])
            print(f'{ odleglosc_2D_3D[ix][0] }|{ odleglosc_2D_3D[ix][1] }')

    else:
        print("\nWybrana opcja nie jest dostępna\n")              



np.savetxt("wspolrzedne.txt", tablica_wynikow, delimiter=',', fmt = ['%12.6f','%14.6f','%15.6f','%15.3f','%15.3f','%15.3f','%15.3f'])