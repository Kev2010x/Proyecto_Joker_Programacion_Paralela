#include <iostream>
#include <thread>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>

void inicializarMatriz(std::vector<float>& matriz, int filas, int columnas) {
    matriz.resize(filas * columnas);
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            matriz[i * columnas + j] = rand() % 100 + 1;
        }
    }
}

void imprimirMatriz(const std::vector<float>& matriz, int filas, int columnas) {
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            std::cout << matriz[i * columnas + j] << " ";
        }
        std::cout << std::endl;
    }
}

void eliminacionAdelante(std::vector<float>& AB, int n, int poscol) {
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

void eliminacionAtras(std::vector<float>& AB, int n, int poscol) {
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

void gaussJordanParallel(std::vector<float>& AB, int n, int numBlocks) {
    std::vector<std::thread> threads;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < numBlocks; ++j) {
            int start = j * (n / numBlocks);
            int end = (j + 1) * (n / numBlocks);

            threads.emplace_back([&AB, n, i, start, end] {
                for (int id = start; id < end; ++id) {
                    eliminacionAdelante(AB, n, i);
                }
            });

            threads.emplace_back([&AB, n, i, start, end] {
                for (int id = start; id < end; ++id) {
                    eliminacionAtras(AB, n, i);
                }
            });
        }

        for (auto &thread : threads) {
            thread.join();
        }

        threads.clear();
    }
}

void gaussJordanSecuencial(std::vector<float>& AB, int n) {
    for (int i = 0; i < n - 1; ++i) {
        for (int id = 0; id < n - 1 - i; ++id) {
            eliminacionAdelante(AB, n, i);
            eliminacionAtras(AB, n, i);
        }
    }
}

int main() {
    const int n = 300;
    srand(static_cast<unsigned>(time(nullptr)));

    std::vector<float> d_AB, d_AB_secuencial;

    inicializarMatriz(d_AB, n + 1, n);
    d_AB_secuencial = d_AB;

    std::cout << "Matriz Inicial (d_AB):" << std::endl;
    imprimirMatriz(d_AB, n + 1, n);

    auto startTimeSecuencial = std::chrono::high_resolution_clock::now();
    gaussJordanSecuencial(d_AB_secuencial, n);
    auto endTimeSecuencial = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sequential_time = endTimeSecuencial - startTimeSecuencial;

    int blocksize = 8;//::min(static_cast<int>(std::thread::hardware_concurrency()), 9);
    int numBlocks = ceil(n / static_cast<float>(blocksize));

    auto startTime = std::chrono::high_resolution_clock::now();
    gaussJordanParallel(d_AB, n, numBlocks);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_time = endTime - startTime;

   // std::cout << "\nMatriz Resultante después de Gauss-Jordan (paralelo):" << std::endl;
   // imprimirMatriz(d_AB, n + 1, n);

   //std::cout << "\nMatriz Resultante después de Gauss-Jordan (secuencial):" << std::endl;
   // imprimirMatriz(d_AB_secuencial, n + 1, n);

    std::cout << "\nSequential Time: " << sequential_time.count() * 1000 << " ms\n";
    std::cout << "Parallel Time: " << parallel_time.count() * 1000 << " ms\n";
    std::cout << "Acceleration: " << sequential_time / parallel_time << "\n";
    std::cout << "Efficiency: " << (sequential_time / parallel_time) / blocksize * 100 << "%\n";

    return 0;
}
