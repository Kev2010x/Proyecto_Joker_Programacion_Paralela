#include <iostream>
#include <cmath>
#include <thread>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>

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
        float piv = AB[pospivot];

        for (int j = pospivot; j < posfinfila; ++j) {
            AB[j] /= piv;
        }

        int posfactor = pospivot + (n + 1) * (id + 1);
        float factor = AB[posfactor];

        for (int j = pospivot; j < posfinfila; ++j) {
            int posactualelim = j + (n + 1) * (id + 1);
            AB[posactualelim] = -1 * factor * AB[j] + AB[posactualelim];
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
        
        AB[posactualelim1] = -1 * factor * AB[pospivot] + AB[posactualelim1];
        AB[posactualelim2] = -1 * factor * AB[pospivot + 1 + poscol] + AB[posactualelim2];
    }
}

void triangSupInf(float *AB, int n, int poscol) {
    for (int idx = 0; idx < n; ++idx) {
        int pospivot = (n + 2) * poscol;
        int posfinfila = (n + 1) * (poscol + 1);
        float piv = AB[pospivot];

        for (int j = pospivot; j < posfinfila; ++j) {
            AB[j] = AB[j] / piv;
        }

        int posfactor = pospivot % (n + 1) + idx * (n + 1);

        if (posfactor != pospivot) {
            float factor = AB[posfactor];

            for (int j = pospivot; j < posfinfila; ++j) {
                int posactualelim = j % (n + 1) + idx * (n + 1);
                AB[posactualelim] = -1 * factor * AB[j] + AB[posactualelim];
            }
        }
    }
}

void resultado(float *AB, int n, float *X) {
    for (int id = 0; id < n; ++id) {
        int posultimacol = (id + 1) * (n + 1) - 1;
        X[id] = AB[posultimacol];
    }
}

void gaussJordan(float *AB, int n, int poscol) {
    eliminacionAdelante(AB, n, poscol);
    eliminacionAtras(AB, n, poscol);
}

int main() {
    const int n = 1000; // Establece el tamaño de la matriz
    srand(static_cast<unsigned>(time(nullptr)));

    float *d_AB = new float[(n + 1) * n];
    float *d_X = new float[n];

    // Inicialización de d_AB
    inicializarMatriz(d_AB, n + 1, n);

    // Imprime la matriz inicial
    std::cout << "Matriz Inicial (d_AB):" << std::endl;
    imprimirMatriz(d_AB, n + 1, n);

    int blocksize = std::thread::hardware_concurrency();
    int numBlocks = ceil(n / static_cast<float>(blocksize));

    auto startTime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n - 1; ++i) {
        std::vector<std::thread> threads;
        for (int j = 0; j < numBlocks; ++j) {
            threads.emplace_back(eliminacionAdelante, d_AB, n, i);
        }

        for (auto &thread : threads) {
            thread.join();
        }
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sequential_time = endTime - startTime;

    startTime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n - 1; ++i) {
        std::vector<std::thread> threads;
        for (int j = 0; j < numBlocks; ++j) {
            threads.emplace_back(eliminacionAtras, d_AB, n, i);  // Corrección: eliminacionAtras en lugar de eliminacionAdelante
        }

        for (auto &thread : threads) {
            thread.join();
        }
    }
    endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_time = endTime - startTime;

    // Imprime información de tiempo
    std::cout << "Tiempo Secuencial Promedio: " << sequential_time.count() * 1000 << " ms\n";
    std::cout << "  Tiempo Paralelo Promedio: " << parallel_time.count() * 1000 << " ms\n";
    std::cout << "Speedup: " << sequential_time / parallel_time << "\n";
    std

::cout << "Eficiencia: " << (sequential_time / parallel_time) / blocksize * 100 << "%\n";

    /*
    // Imprime el vector de soluciones
    std::cout << "Vector de Soluciones (X): ";
    for (int i = 0; i < n; ++i) {
        std::cout << d_X[i] << " ";
    }
    std::cout << std::endl;
    */

    // Libera la memoria y utiliza los resultados en d_X
    delete[] d_AB;
    delete[] d_X;

    return 0;
}

