
#include <iostream>
#include <cmath>
#include <thread>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <future>
#include <mutex>
#include <condition_variable>

// Mutex y variable de condición para la sincronización
std::mutex mutexAB;
std::condition_variable cv;

void inicializarMatriz(float *matriz, int filas, int columnas) {
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            matriz[i * columnas + j] = rand() % 100 + 1;  // Valor aleatorio entre 1 y 100
        }
    }
}

void imprimirMatriz(const float *matriz, int filas, int columnas) {
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
        float piv;

        // Sección crítica: leer y actualizar AB
        {
            std::unique_lock<std::mutex> lock(mutexAB);
            piv = AB[pospivot];

            for (int j = pospivot; j < posfinfila; ++j) {
                AB[j] /= piv;
            }
        }  // Fin de la sección crítica

        int posfactor = pospivot + (n + 1) * (id + 1);
        float factor;

        // Sección crítica: leer AB
        {
            std::unique_lock<std::mutex> lock(mutexAB);
            factor = AB[posfactor];
        }  // Fin de la sección crítica

        // Sección crítica: actualizar AB
        {
            std::unique_lock<std::mutex> lock(mutexAB);
            for (int j = pospivot; j < posfinfila; ++j) {
                int posactualelim = j + (n + 1) * (id + 1);
                AB[posactualelim] = -1 * factor * AB[j] + AB[posactualelim];
            }
        }  // Fin de la sección crítica
    }
}

void eliminacionAtras(float *AB, int n, int poscol) {
    for (int id = 0; id < n - 1 - poscol; ++id) {
        int pospivot = (n + 2) * (n - 1 - poscol);

        if (poscol == 0) {
            float pivot;
            
            // Sección crítica: leer y actualizar AB
            {
                std::unique_lock<std::mutex> lock(mutexAB);
                pivot = AB[pospivot];
                AB[pospivot] = AB[pospivot] / pivot;
                AB[pospivot + 1] = AB[pospivot + 1] / pivot;
            }  // Fin de la sección crítica
        }

        float factor;
        int posactualelim1 = pospivot - (n + 1) * (id + 1);
        int posactualelim2 = pospivot - (n + 1) * (id + 1) + 1 + poscol;

        // Sección crítica: leer AB
        {
            std::unique_lock<std::mutex> lock(mutexAB);
            factor = AB[pospivot - (n + 1) * (id + 1)];
        }  // Fin de la sección crítica
        
        // Sección crítica: actualizar AB
        {
            std::unique_lock<std::mutex> lock(mutexAB);
            AB[posactualelim1] = -1 * factor * AB[pospivot] + AB[posactualelim1];
            AB[posactualelim2] = -1 * factor * AB[pospivot + 1 + poscol] + AB[posactualelim2];
        }  // Fin de la sección crítica
    }
}

void gaussJordanParallel(float *AB, int n, int numBlocks) {
    std::vector<std::future<void>> futures;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < numBlocks; ++j) {
            int start = j * (n / numBlocks);
            int end = (j + 1) * (n / numBlocks);

            // Lanzar hilos para la eliminación hacia adelante
            futures.emplace_back(std::async(std::launch::async, [AB, n, i, start, end] {
                for (int id = start; id < end; ++id) {
                    eliminacionAdelante(AB, n, i);
                }
            }));

            // Lanzar hilos para la eliminación hacia atrás
            futures.emplace_back(std::async(std::launch::async, [AB, n, i, start, end] {
                for (int id = start; id < end; ++id) {
                    eliminacionAtras(AB, n, i);
                }
            }));
        }

        // Esperar a que todos los bloques lanzados en esta iteración completen su trabajo.
        for (auto &future : futures) {
            future.wait();
        }

        futures.clear();
    }
}

// Función para realizar el método de eliminación de Gauss-Jordan de manera secuencial
void gaussJordanSecuencial(float *AB, int n) {
    for (int i = 0; i < n - 1; ++i) {
        for (int id = 0; id < n - 1 - i; ++id) {
            eliminacionAdelante(AB, n, i);
            eliminacionAtras(AB, n, i);
        }
    }
}

int main() {
    const int n = 800; // Establece el tamaño de la matriz
    srand(static_cast<unsigned>(time(nullptr)));

    float *d_AB = new float[(n + 1) * n];
    float *d_AB_secuencial = new float[(n + 1) * n]; // Copia para el método secuencial

    // Inicialización de d_AB
    inicializarMatriz(d_AB, n + 1, n);
    memcpy(d_AB_secuencial, d_AB, sizeof(float) * (n + 1) * n);  // Copia para el método secuencial

    // Imprime la matriz inicial
    std::cout << "Matriz Inicial (d_AB):" << std::endl;
    imprimirMatriz(d_AB, n + 1, n);

    // Medir tiempo secuencial
    auto startTimeSecuencial = std::chrono::high_resolution_clock::now();
    gaussJordanSecuencial(d_AB_secuencial, n);
    auto endTimeSecuencial = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sequential_time = endTimeSecuencial - startTimeSecuencial;
   
    // medir tiempo de ejecución paralella
    int blocksize = std::thread::hardware_concurrency();
    int numBlocks = ceil(n / static_cast<float>(blocksize));

    
    auto startTime = std::chrono::high_resolution_clock::now();
    gaussJordanParallel(d_AB, n, numBlocks);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_time = endTime - startTime;

    // Imprime información de tiempo
   // std::cout << "Tiempo Paralelo Promedio: " << parallel_time.count() * 1000 << " ms\n";

    // Imprime información de tiempo
    std::cout << "Tiempo Secuencial: " << sequential_time.count() * 1000 << " ms\n";
    std::cout << "Tiempo Paralelo: " << parallel_time.count() * 1000 << " ms\n";
    std::cout << "Aceleración: " << sequential_time / parallel_time << "\n";
    std::cout << "Eficiencia: " << (sequential_time / parallel_time) / blocksize * 100 << "%\n";

    // Libera la memoria
    delete[] d_AB;
    delete[] d_AB_secuencial;

    return 0;
}


