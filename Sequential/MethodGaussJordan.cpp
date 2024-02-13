#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <random>

using namespace std;

void gauss_jordan(double M[500][500], int n){
    double mayor;
    int indice;
    double aux;
    double pivote;

    for(int k=0;k<n;k++){
        mayor=abs(M[k][k]);
        indice=k;
        for(int l=k+1;l<n;l++){
            if(mayor<abs(M[l][k]))
            {
                mayor=abs(M[l][k]);
                indice=l;
            }
        }
        if(k!=indice){
            for(int i =0;i<n+1;i++){
                aux=M[k][i];
                M[k][i]=M[indice][i];
                M[indice][i]=aux;
            }
        }
        if(M[k][k]==0){
            cout<<"El sistema no tiene solución"<<endl;
            exit(0);
        }
        else{
            for(int i=0;i<n;i++){
                if(i!=k){
                    pivote=-M[i][k];
                    for(int j=k;j<n+1;j++){
                        M[i][j]=M[i][j]+pivote*M[k][j]/M[k][k];
                    }
                }
                else{
                    pivote=M[k][k];
                    for(int j=k;j<n+1;j++){
                        M[i][j]=M[i][j]/pivote;
                    }
                }
            }
        }
    }
}

void ingresar_coeficientes_aleatorios(double M[500][500], int n){
    // Generador de números aleatorios
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(-10.0, 10.0); // Números aleatorios en el rango [-10.0, 10.0]

    for(int i=0; i<n; i++){
        for(int j=0; j<n+1; j++){
            M[i][j] = dis(gen); // Llenar la matriz con números aleatorios
        }
    }
}

void mostrar_matriz(double M[500][500], int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n+1;j++){
            cout<<M[i][j]<<" ";
        }
        cout<<endl;
    }
}

void mostrar_matriz_aumentada(double M[500][500], int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n+1;j++){
            cout<<M[i][j]<<" ";
        }
        cout<<endl;
    }
}

void mostrar_resultados(double M[500][500], int n){
    cout<<"Los resultados son: "<<endl;
    for(int i=0;i<n;i++){
        cout<<"x"<<i+1<<" = "<<M[i][n]<<endl;
    }
}

int main(){
    double M[500][500];
    int n;
    cout << "Ingrese el número de ecuaciones: ";
    cin >> n;

    ingresar_coeficientes_aleatorios(M, n);

    cout << "| Matriz |" << endl;
    mostrar_matriz(M, n);

    auto start_time = chrono::high_resolution_clock::now();

    gauss_jordan(M, n);

    auto end_time = chrono::high_resolution_clock::now();

    cout << "| Matriz reducida |" << endl;
    mostrar_matriz_aumentada(M, n);

    cout << "| Resultados |" << endl;
    mostrar_resultados(M, n);

    auto duration_ms = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cout << "Tiempo de ejecución: " << duration_ms.count() << " ms" << endl;

    return 0;
}
