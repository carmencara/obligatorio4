//---------------------------------------------------------------------
//                                                                    |
//         ALGORITMO PARA ENVIAR UN COHETE A LA LUNA USANDO RK4       |
//                              Objetivos:                            |
//                                                                    |
//   1. Resolver las ecuaciones de movimiento para enviar el cohete   |
//      a la Luna.                                                    |
//   2. Comprobar que se conserva H'=H-w*p_phi.                       | 
//                                                                    |
//---------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <fstream>      // Para trabajar con ficheros

#define PI 3.14159265
#define G 6.67E-11      // Constante de gravitación universal (en N·m^2/kg^2)
#define M_T 5.9736E24   // Masa de la Tierra (en kg)
#define M_L 0.07349E24  // Masa de la Luna (en kg)
#define d_TL 3.844E8    // Distancia fija entre la Tierra y la Luna (en m)
#define w 2.6617E-6     // Velocidad angular constante de la Luna (en s^-1)
#define R_T 6.37816E6   // Radio de la Tierra (en m)
#define R_L 1.7374E6    // Radio de la Luna (en m)

using namespace std;

// FUNCIONES QUE SE VAN A UTILIZAR


/*---------------------------------------------------------------------
|                           FUNCIÓN PRINCIPAL                         |
---------------------------------------------------------------------*/
int main()
{
    // ------------------ DECLARACIÓN DE VARIABLES --------------------
    int i,j;                // Contadores
    
    double h;               // Paso
    double delta;           // Delta=G*M_T/(d_TL)^3
    double mu;              // mu=M_L/M_T
    double rtilde_prima;    // r~'=sqrt(1+r^2-2*r*cos(phi-wt))
    double x[4];            // Vector donde se guardan: r~, phi, p~_r, p~_phi
                            // x[0]=r~=r/d_TL ; x[1]=phi ; x[2]=p~_r ; x[3]=p~_phi
    double v_esc;           // Velocidad de escape de la Tierra

    ofstream fich_posicion; // Fichero para guardar la posición de los cuerpos

    // Abrir el fichero de la posición de los cuerpos en cada instante
    fich_posicion.open("datos.txt");

    // ------------------------ INICIALIZACIÓN ------------------------
    h = 1.0;
    delta = G*M_T/pow(d_TL,3);
    mu = M_L / M_T;
    v_esc = // 11.2 km/s

    // Inicializar las componentes de x[4]
    x[0] = R_T/d_TL;        // r~(0)=R_T/d_TL El cohete sale de la superficie de la Tierra
    x[1] = PI/2;            // El cohete sale a lo largo del eje y (phi=90º)
    x[2] = 

    fich_posicion.close();
    return 0;
}