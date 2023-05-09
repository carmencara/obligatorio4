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
#define R_T 6.378160E6  // Radio de la Tierra (en m)
#define R_L 1.7374E6    // Radio de la Luna (en m)

using namespace std;

// FUNCIONES QUE SE VAN A UTILIZAR
double RK4 (double (&x)[4], double (&k)[4][4], double h, double delta, double mu, double r_prima, double t);
double Derivada_r (double x0, double x1, double x2, double x3);
double Derivada_phi (double x0, double x1, double x2, double x3);
double Derivada_pr (double x0, double x1, double x2, double x3, double r_prima, double delta, double mu, double t);
double Derivada_pphi (double x0, double x1, double x2, double x3, double r_prima, double delta, double mu, double t);


/*---------------------------------------------------------------------
|                           FUNCIÓN PRINCIPAL                         |
---------------------------------------------------------------------*/
int main()
{
    // ------------------ DECLARACIÓN DE VARIABLES --------------------
    int i,j;                // Contadores
    int iteraciones;        // Número de iteraciones del algoritmo RK4
    
    double h;               // Paso
    double t;               // Tiempo
    double delta;           // Delta=G*M_T/(d_TL)^3
    double mu;              // mu=M_L/M_T
    double r_prima;         // r~'=sqrt(1+r^2-2*r*cos(phi-wt))
    double x[4];            // Vector: x[0]=r~=r/d_TL ; x[1]=ϕ ; x[2]=p~_r ; x[3]=p~_ϕ
    double k[4][4];         // Matriz 4x4 con las k para las 4 ecuaciones diferenciales
    double r0;              // r~ inicial: el cohete sale de la superficie de la Tierra
    double v_esc;           // Velocidad de escape de la Tierra (11.19 km/s)
    double theta0;          // Ángulo inicial entre el eje x y la velocidad de lanzamiento
    double phi0;            // Ángulo inicial entre el eje x y el punto de lanzamiento

    ofstream fich_posicion; // Fichero para guardar la posición de los cuerpos

    // Abrir el fichero de la posición de los cuerpos en cada instante
    fich_posicion.open("datos.txt");

    // ------------------------ INICIALIZACIÓN ------------------------
    h = 1.0;
    t = 0.0;
    delta = G*M_T/pow(d_TL,3);
    mu = M_L / M_T;

    r0 = R_T/d_TL;
    v_esc = 11190/d_TL;
    theta0 = 0.45;
    phi0 = 0.45;

    // Inicializar las componentes de x[4]
    x[0] = r0;                        // r~(0)=R_T/d_TL
    x[1] = phi0;                      // El cohete sale a lo largo del eje y (ϕ=90º)
    x[2] = v_esc*cos(theta0-phi0);    // p~_r(0)=v_esc*cos(theta0-phi0)
    x[3] = r0*v_esc*sin(theta0-phi0); // p~_ϕ(0)=r~(0)*v_esc*sin(theta0-phi0)

    // Guardar en el fichero la situación inicial de la Luna y el cohete (la Tierra permanece en (0,0))
    fich_posicion << 0 << ", " << 0 << endl; // Tierra
    fich_posicion << cos(w*t) << ", " << sin(w*t) << endl; // Luna
    fich_posicion << x[0]*cos(x[1]) << ", " << x[0]*sin(x[1]) << endl << endl; // Cohete


    // ------------------------- AGORITMO RK4 -------------------------
    iteraciones = 200000;
    for(i=0; i<=iteraciones; i++)
    {
        // 1. Calcular r~'=sqrt(1+r~^2-2*r~*cos(ϕ-wt))
        r_prima = sqrt(1+x[0]*x[0]-2*x[0]*cos(x[1]-w*t));
        // 2. Aplicar el algoritmo de Runge-Kutta de orden 4 y actualizar el tiempo t=t+h
        t = t + RK4(x, k, h, delta, mu, r_prima, t);
        // 3. Escribir las posiciones de la Tierra, Luna y cohete (no en todas las iteraciones)
        if(i%500==0)
        {
            fich_posicion << 0 << ", " << 0 << endl; // Tierra
            fich_posicion << cos(w*t) << ", " << sin(w*t) << endl; // Luna
            fich_posicion << x[0]*cos(x[1]) << ", " << x[0]*sin(x[1]) << endl << endl; // Cohete
        }
    }

    fich_posicion.close();
    return 0;
}


/*---------------------------------------------------------------------
|                               ALGORITMO                             |
---------------------------------------------------------------------*/

// FUNCIÓN RK4: aplica el algoritmo de Runge-Kutta de orden 4 (calcula k[][] y actualiza x[])
// y devuelve el paso
double RK4 (double (&x)[4], double (&k)[4][4], double h, double delta, double mu, double r_prima, double t)
{   
    int j; // Contador

    // Calcular las componentes de k[][]
    k[0][0]=h*Derivada_r(x[0], x[1], x[2], x[3]);
    k[0][1]=h*Derivada_phi(x[0], x[1], x[2], x[3]);
    k[0][2]=h*Derivada_pr(x[0], x[1], x[2], x[3], r_prima, delta, mu, t);
    k[0][3]=h*Derivada_pphi(x[0], x[1], x[2], x[3], r_prima, delta, mu, t);

    k[1][0]=h*Derivada_r(x[0]+k[0][0]/2, x[1]+k[0][1]/2, x[2]+k[0][2]/2, x[3]+k[0][3]/2);
    k[1][1]=h*Derivada_phi(x[0]+k[0][0]/2, x[1]+k[0][1]/2, x[2]+k[0][2]/2, x[3]+k[0][3]/2);
    k[1][2]=h*Derivada_pr(x[0]+k[0][0]/2, x[1]+k[0][1]/2, x[2]+k[0][2]/2, x[3]+k[0][3]/2, r_prima, delta, mu, t+h/2);
    k[1][3]=h*Derivada_pphi(x[0]+k[0][0]/2, x[1]+k[0][1]/2, x[2]+k[0][2]/2, x[3]+k[0][3]/2, r_prima, delta, mu, t+h/2);     
        
    k[2][0]=h*Derivada_r(x[0]+k[1][0]/2, x[1]+k[1][1]/2, x[2]+k[1][2]/2, x[3]+k[1][3]/2);
    k[2][1]=h*Derivada_phi(x[0]+k[1][0]/2, x[1]+k[1][1]/2, x[2]+k[1][2]/2, x[3]+k[1][3]/2);
    k[2][2]=h*Derivada_pr(x[0]+k[1][0]/2, x[1]+k[1][1]/2, x[2]+k[1][2]/2, x[3]+k[1][3]/2, r_prima, delta, mu, t+h/2);
    k[2][3]=h*Derivada_pphi(x[0]+k[1][0]/2, x[1]+k[1][1]/2, x[2]+k[1][2]/2, x[3]+k[1][3]/2, r_prima, delta, mu, t+h/2);

    k[3][0]=h*Derivada_r(x[0]+k[2][0], x[1]+k[2][1], x[2]+k[2][2], x[3]+k[2][3]);
    k[3][1]=h*Derivada_phi(x[0]+k[2][0], x[1]+k[2][1], x[2]+k[2][2], x[3]+k[2][3]);
    k[3][2]=h*Derivada_pr(x[0]+k[2][0], x[1]+k[2][1], x[2]+k[2][2], x[3]+k[2][3], r_prima, delta, mu, t+h);
    k[3][3]=h*Derivada_pphi(x[0]+k[2][0], x[1]+k[2][1], x[2]+k[2][2], x[3]+k[2][3], r_prima, delta, mu, t+h); 

    // Actualizar las componentes de x[]
    for(j=0; j<4; j++) x[j] = x[j] + 1./6*(k[0][j]+2*k[1][j]+2*k[2][j]+k[3][j]);

    return h;
}

/*---------------------------------------------------------------------
|                       ECUACIONES DEL MOVIMIENTO                     |
---------------------------------------------------------------------*/
// Se les pasa como argumento las componentes del vector x[] necesarias:
// x[0]=x0=r~, x[1]=x1=ϕ, x[2]=x2=p~_r, x[3]=x3=p~_ϕ y otras variables

// FUNCIÓN Derivada_r: devuelve la primera ecuación de movimiento: r~^punto=p~_r
double Derivada_r (double x0, double x1, double x2, double x3)
{
    return x2;
}

// FUNCIÓN Derivada_phi: devuelve la segunda ecuación de movimiento: ϕ^punto=p~_ϕ/r~^2
double Derivada_phi (double x0, double x1, double x2, double x3)
{
    return x3/(x0*x0);
}

// FUNCIÓN Derivada_pr: devuelve la tercera ecuación de movimiento: p~_r^punto=...
double Derivada_pr (double x0, double x1, double x2, double x3, double r_prima, double delta, double mu, double t)
{
    return x3*x3/pow(x0,3) - delta*((1/(x0*x0))+(mu*(x0-cos(x1-w*t))/pow(r_prima,3)));
}

// FUNCIÓN Derivada_phi: devuelve la tercera ecuación de movimiento: p~_ϕ^punto=...
double Derivada_pphi (double x0, double x1, double x2, double x3, double r_prima, double delta, double mu, double t)
{
    return -1.0*delta*mu*x0*sin(x1-w*t)/pow(r_prima,3);
}