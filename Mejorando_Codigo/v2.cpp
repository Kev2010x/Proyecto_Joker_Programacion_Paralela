// en proceso de mejora
#include <iostream>
#include <cmath>
#include <thread>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <mutex>
#include <omp.h>
#include <condition_variable>

// Mutex para garantizar acceso seguro a la impresión en consola
std::mutex printMutex;

// Barrera para sincronizar los hilos
std::mutex barrierMutex;
std::condition_variable barrierCV;
int barrierCount = 0;
const int numThreads = std::thread::hardware_concurrency();

void inicializarMatriz(float *matriz, int filas, int columnas) {
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            matriz[i * columnas + j] = rand() % 100 + 1;  // Valor aleatorio entre 1 y 100
        }
    }
}

void imprimirMatriz(const float *matriz, int filas, int columnas) {
    std::lock_guard<std::mutex> lock(printMutex);  // Mutex para impresión en consola
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            std::cout << matriz[i * columnas + j] << " ";
        }
        std::cout << std::endl;
    }
}

void eliminacionAdelante(float *AB, int n, int poscol) {
    for (int id = 0; id < n - 1 - poscol; ++id) {
        int pospivot = (n + 2) * poscol;
        int posfinfila = (n + 1) * (poscol + 1);
        float piv = AB[pospivot];

        for (int j = pospivot; j < posfinfila; ++j) {
            AB[j] /= piv;
        }

        #pragma omp parallel for num_threads(numThreads)
        for (int j = pospivot + 1; j < posfinfila; ++j) {
            int posactualelim = j + (n + 1) * (id + 1);
            AB[posactualelim] -= AB[j];
        }
    }
}

void eliminacionAtras(float *AB, int n, int poscol) {
    for (int id = 0; id < n - 1 - poscol; ++id) {
        int pospivot = (n + 2) * (n - 1 - poscol);

        if (poscol == 0) {
            float pivot = AB[pospivot];
            AB[pospivot] = AB[pospivot] / pivot;
            AB[pospivot + 1] = AB[pospivot + 1] / pivot;
        }

        float factor = AB[pospivot - (n + 1) * (id + 1)];
        int posactualelim1 = pospivot - (n + 1) * (id + 1);
        int posactualelim2 = pospivot - (n + 1) * (id + 1) + 1 + poscol;

        AB[posactualelim1] -= factor * AB[pospivot];
        AB[posactualelim2] -= factor * AB[pospivot + 1 + poscol];
    }
}

void gaussJordan(float *AB, int n, int poscol) {
    eliminacionAdelante(AB, n, poscol);
    eliminacionAtras(AB, n, poscol);
}

void barrier(int id) {
    std::unique_lock<std::mutex> lock(barrierMutex);
    barrierCount++;
    if (barrierCount == numThreads) {
        barrierCount = 0;
        barrierCV.notify_all();
    } else {
        barrierCV.wait(lock, [&] { return barrierCount == 0; });
    }
}

int main() {
    const int n = 20; // Establece el tamaño de la matriz
    srand(static_cast<unsigned>(time(nullptr)));

    float *d_AB = new float[(n + 1) * n];

    // Inicialización de d_AB
    inicializarMatriz(d_AB, n + 1, n);

    // Imprimir la matriz inicial
    imprimirMatriz(d_AB, n + 1, n);

    auto startTime = std::chrono::high_resolution_clock::now();

    // Proceso de eliminación de Gauss-Jordan en paralelo
    #pragma omp parallel num_threads(numThreads)
    {
        for (int i = 0; i < n - 1; ++i) {
            gaussJordan(d_AB, n, i);
            barrier(i);
        }
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_time = endTime - startTime;

    // Imprimir información de tiempo
    std::cout << "Tiempo Paralelo: " << parallel_time.count() * 1000 << " ms\n";

    // Imprimir el vector solución
    std::cout << "Vector Solución (X): ";
    for (int i = 0; i < n; ++i) {
        std::cout << d_AB[i + (n + 1) * (n - 1)] << " ";
    }
    std::cout << std::endl;

    // Liberar memoria
    delete[] d_AB;

    return 0;
}
