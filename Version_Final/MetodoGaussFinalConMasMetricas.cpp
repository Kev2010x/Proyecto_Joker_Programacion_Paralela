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
#include <atomic>

std::mutex mutexAB;
std::condition_variable cv;

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

void eliminacionAdelante(std::vector<float>& AB, int n, int poscol, std::atomic<int>& total_eliminations) {
    for (int id = 0; id < n - 1 - poscol; ++id) {
        int pospivot = (n + 2) * poscol;
        int posfinfila = (n + 1) * (poscol + 1);
        float piv;

        {
            std::unique_lock<std::mutex> lock(mutexAB);
            piv = AB[pospivot];

            for (int j = pospivot; j < posfinfila; ++j) {
                AB[j] /= piv;
            }
        }

        int posfactor = pospivot + (n + 1) * (id + 1);
        float factor;

        {
            std::unique_lock<std::mutex> lock(mutexAB);
            factor = AB[posfactor];
        }

        {
            std::unique_lock<std::mutex> lock(mutexAB);
            for (int j = pospivot; j < posfinfila; ++j) {
                int posactualelim = j + (n + 1) * (id + 1);
                AB[posactualelim] = -1 * factor * AB[j] + AB[posactualelim];
            }
        }
        total_eliminations += (n - poscol);
    }
}

void eliminacionAtras(std::vector<float>& AB, int n, int poscol, std::atomic<int>& total_eliminations) {
    for (int id = 0; id < n - 1 - poscol; ++id) {
        int pospivot = (n + 2) * (n - 1 - poscol);

        if (poscol == 0) {
            float pivot;

            {
                std::unique_lock<std::mutex> lock(mutexAB);
                pivot = AB[pospivot];
                AB[pospivot] = AB[pospivot] / pivot;
                AB[pospivot + 1] = AB[pospivot + 1] / pivot;
            }
        }

        float factor;
        int posactualelim1 = pospivot - (n + 1) * (id + 1);
        int posactualelim2 = pospivot - (n + 1) * (id + 1) + 1 + poscol;

        {
            std::unique_lock<std::mutex> lock(mutexAB);
            factor = AB[pospivot - (n + 1) * (id + 1)];
        }

        {
            std::unique_lock<std::mutex> lock(mutexAB);
            AB[posactualelim1] = -1 * factor * AB[pospivot] + AB[posactualelim1];
            AB[posactualelim2] = -1 * factor * AB[pospivot + 1 + poscol] + AB[posactualelim2];
        }
        total_eliminations += (n - poscol);
    }
}

void gaussJordanParallel(std::vector<float>& AB, int n, int numBlocks, std::atomic<int>& total_eliminations) {
    std::vector<std::future<void>> futures;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < numBlocks; ++j) {
            int start = j * (n / numBlocks);
            int end = (j + 1) * (n / numBlocks);

            futures.emplace_back(std::async(std::launch::async, [&AB, n, i, start, end, &total_eliminations] {
                for (int id = start; id < end; ++id) {
                    eliminacionAdelante(AB, n, i, total_eliminations);
                }
            }));

            futures.emplace_back(std::async(std::launch::async, [&AB, n, i, start, end, &total_eliminations] {
                for (int id = start; id < end; ++id) {
                    eliminacionAtras(AB, n, i, total_eliminations);
                }
            }));
        }
    }

    for (auto &future : futures) {
        future.wait();
    }
    futures.clear();  // Mueve esto fuera del bucle externo
}

void gaussJordanSecuencial(std::vector<float>& AB, int n, std::atomic<int>& total_eliminations) {
    for (int i = 0; i < n - 1; ++i) {
        for (int id = 0; id < n - 1 - i; ++id) {
            eliminacionAdelante(AB, n, i, total_eliminations);
            eliminacionAtras(AB, n, i, total_eliminations);
        }
    }
}

int main() {
    const int n = 100;
    srand(static_cast<unsigned>(time(nullptr)));

    std::vector<float> d_AB, d_AB_secuencial;

    inicializarMatriz(d_AB, n + 1, n);
    d_AB_secuencial = d_AB;

    std::cout << "Matriz Inicial (d_AB):" << std::endl;
    imprimirMatriz(d_AB, n + 1, n);

    std::atomic<int> total_eliminations_secuencial(0);
    std::atomic<int> total_eliminations_paralelo(0);

    auto startTimeSecuencial = std::chrono::high_resolution_clock::now();
    gaussJordanSecuencial(d_AB_secuencial, n, total_eliminations_secuencial);
    auto endTimeSecuencial = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sequential_time = endTimeSecuencial - startTimeSecuencial;

    int blocksize = std::min(static_cast<int>(std::thread::hardware_concurrency()), 12);
    int numBlocks = ceil(n / static_cast<float>(blocksize));

    std::cout << "Number of threads used: " << blocksize << std::endl;

    auto startTime = std::chrono::high_resolution_clock::now();
    gaussJordanParallel(d_AB, n, numBlocks, total_eliminations_paralelo);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_time = endTime - startTime;

    std::cout << "Sequential Time: " << sequential_time.count() * 1000 << " ms\n";
    std::cout << "Parallel Time: " << parallel_time.count() * 1000 << " ms\n";
    std::cout << "Acceleration: " << sequential_time / parallel_time << "\n";
    std::cout << "Efficiency: " << (sequential_time / parallel_time) / blocksize * 100 << "%\n";

    double throughput_secuencial = total_eliminations_secuencial / sequential_time.count();
    double throughput_paralelo = total_eliminations_paralelo / parallel_time.count();
    double latency = 1 / throughput_paralelo;
    std::cout << "Sequential Performance: " << throughput_secuencial << " Operations/Second\n";
    std::cout << "Parallel Performance: " << throughput_paralelo << " Operations/Second\n";
    std::cout << "Latency: " <<  latency  << " Operations/Second\n";
    
    return 0;
}
