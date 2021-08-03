/// Autor: Tomasz Boron

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include<fstream>
#include<stdlib.h>
#include<unistd.h>

using namespace std;
using namespace Eigen;

/// omega reprezentowana przez przedzial [L,R]
/// H wyliczane na start
/// N to ilosc przedzialow
double L=0;
double R=2;
double N=11;
double H=(R-L)/N;

double ei(double elem,int i){
    double xi=L+i*H;
    if(elem<=L || elem>R) return 0;
    if(elem<=xi && elem>=xi-H) return (elem-(xi-H))/H;
    if(elem>=xi && elem<=xi+H) return (xi+H-elem)/H;
    return 0;
}

double eiPrim(double elem, int i){ /// poniewaz funkcje liniowe to wyznaczamy wspólczynnik przy x jako pochodna
    double xi=L+i*H;
    if(elem<=L || elem>R) return 0;
    if(elem<=xi && elem>=xi-H) return 1/H;
    if(elem>=xi && elem<=xi+H) return -1/H;
    return 0;
}

double sinus(double x, int anything){ /// wlasny, aby argumenty sie spasowaly, drugi argumeny nic nie wnosi
    return sin(x);
}

double integral(double (*u)(double,int),double (*v)(double,int),double a,double b,int i, int j){ /// calkowanie iloczynu funkcji metoda trapezowa
    double score=0;
    int n=10000;
    double h=(b-a)/double(n);
    for(int k=1;k<=n;k++){
        double p1=(*u)(a+(k-1)*h,i)*(*v)(a+(k-1)*h,j);
        double p2=(*u)(a+k*h,i)*(*v)(a+k*h,j);
        score+=(p1+p2)*h/2;
    }
    return score;
}

double Buv(double (*u)(double,int),double (*v)(double,int),double (*u_prim)(double,int),double (*v_prim)(double,int),double left,double right,int i,int j){
    return -(*v)(right,j)*(*u)(right,i)+integral(u_prim,v_prim,left,right,i,j)-integral(u,v,left,right,i,j);
}

double Lv(double (*v)(double,int),double left,double right,int i){
    return integral(v,sinus,left,right,i,i);
}

double result(double elem,VectorXf u){
    double score=0;
    for(int i=0;i<=N;i++){
        score+=u(i)*ei(elem,i);
    }
    return score;
}

int main()
{
    MatrixXf A((int)N+1,(int)N+1);
    for(int i=1;i<=N;i++){      /// wymuszenie warunku Dirichleta (jak w filmiku lab9) cz I
        A(0,i)=0;
        A(i,0)=0;
    }
    A(0,0)=1;
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            if(abs(i-j)<=1){    /// poniewaz dalej niz po sasiedzku te piramidki niezerowe beda rozlaczne, iloczyn pod calkami bedzie wszedzie 0
                A(j,i)=Buv(ei,ei,eiPrim,eiPrim,L,R,i,j);
            }
            else A(j,i)=0;
        }
    }
    cout << "Macierz A:\n" << A << endl;
    VectorXf b((int)N+1);
    b(0)=0;                     /// wymuszenie warunku Dirichleta (jak w filmiku lab 9) cz II
    for(int i=1;i<=N;i++){
        b(i)=Lv(ei,L,R,i);
    }
    cout << "Wektor b:\n" << b << endl;

    /// colPivHouseholder - dekompozycja macierzy
    /// rozklada macierz na macierze P,Q,R takie ze
    /// P - macierz permutacji (w kazdym wierszu i kolumnie dokladnie jedna 1,reszta 0)
    /// Q - macierz unitarna
    /// R - macierz trojkatna gorna
    VectorXf x=A.colPivHouseholderQr().solve(b);
    cout << "Rozwiazaniem jest wektor x:\n";
    for(int i=0;i<=N;i++) cout << i << ". " << x(i) << endl;
    cout << "Wcisnij enter aby wyjsc\n";
    FILE *pipe = popen("gnuplot", "w");
    if (pipe != NULL)
    {
        fprintf(pipe, "set xlabel 'os X'\n");
        fprintf(pipe, "set ylabel 'os Y'\n");
        fprintf(pipe, "set grid linestyle 3\n");

        fprintf(pipe, "plot '-' with linespoints linetype 1 linewidth 1\n");

        for(double i = L; i < R; i+=0.001) fprintf(pipe, "%f %f\n",i, result(i,x));
        fprintf(pipe, "%s\n", "e");
        fflush(pipe);
        cin.clear();
        cin.get();
        pclose(pipe);
    }
    else cout << "Blad otwarcia\n" << endl;
    /// biblioteka Eigen: https://eigen.tuxfamily.org/
    /// gnuplot:
    ///     https://stackoverflow.com/questions/9349538/gnuplot-c-from-windows-command-window-opens-and-closes
    ///     https://cpp0x.pl/artykuly/Inne-artykuly/C++-Programistyczne-ustawianie-wykresow-dla-programu-Gnuplot/52
    /// przekazanie funkcji do innej funkcji: https://www.quantstart.com/articles/Function-Objects-Functors-in-C-Part-1/
}
