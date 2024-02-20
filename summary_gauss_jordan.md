# Gauss Jordan Méthod

`Librerias`

#include <thread: Permite trabajar con concurrencia y multihilo en C++.

#include <vector: Proporciona la clase vector para trabajar con arreglos dinámicos en C++.

#include <cstdlib: Proporciona funciones relacionadas con el sistema, como la generación de números aleatorios. 

#include <ctime: Ofrece funciones para trabajar con el tiempo en C++.

#include <chrono: Proporciona clases y funciones para medir el tiempo en C++. Se utiliza para calcular el tiempo de ejecución de ciertas operaciones.

#include <future: Permite trabajar con resultados futuros de operaciones asíncronas en C++. Se utiliza para realizar operaciones asíncronas en programación paralela.

#include <mutex: Proporciona clases y funciones para trabajar con mutex en C++. Se utiliza para sincronizar el acceso concurrente a secciones críticas mediante bloqueos de mutex.

#include <condition_variable: Proporciona clases y funciones para trabajar con variables de condición en C++.
```c++
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

```
 Mutex y condición variable para garantizar la sincronización en el acceso a la matriz AB
```c++
    std::mutex mutexAB;
    std::condition_variable cv;
```

inicializarMatriz: Esta función se encarga de llenar una matriz cuadrada de tamaño (filas + 1) x columnas con valores aleatorios entre 1 y 100. La matriz representa el sistema de ecuaciones aumentado.
```c++
void inicializarMatriz(std::vector<float>& matriz, int filas, int columnas) {
    matriz.resize(filas * columnas);
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            matriz[i * columnas + j] = rand() % 100 + 1;
        }
    }
}
```



Impresión de la matriz
```c++
void imprimirMatriz(const std::vector<float>& matriz, int filas, int columnas) {
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            std::cout << matriz[i * columnas + j] << " ";
        }
        std::cout << std::endl;
    }
}
```
## Funcion EliminacionAdelante

La función `eliminacionAdelante` realiza la eliminación hacia adelante en el método de Gauss-Jordan. Se utiliza un mutex para garantizar el acceso exclusivo a la sección crítica de la matriz aumentada AB, evitando así condiciones de carrera en operaciones concurrentes. Itera sobre las filas y normaliza la fila actual dividiendo por el pivote. Luego, resta múltiplos de esta fila a las filas siguientes.  //  en una columna específica
```c++
void eliminacionAdelante(std::vector<float>& AB, int n, int poscol) {
    for (int id = 0; id < n - 1 - poscol; ++id) { // calculo de indices
        int pospivot = (n + 2) * poscol; // Índice de la posición del pivote en la fila actual.
        int posfinfila = (n + 1) * (poscol + 1); // Índice que marca el final de la fila actual.
        float piv;

        {
            // Bloquear el mutex para garantizar acceso exclusivo a la sección crítica
            std::unique_lock<std::mutex> lock(mutexAB);
            piv = AB[pospivot];

            // Realizar operaciones en la matriz dentro de la sección crítica
            for (int j = pospivot; j < posfinfila; ++j) {
                AB[j] /= piv; // division de toda la fila/pivote
            }
        }
        // Índice de la posición del factor en la fila actual.
        int posfactor = pospivot + (n + 1) * (id + 1);
        float factor;

        {
            // Bloquear el mutex para garantizar acceso exclusivo a la sección crítica

            //Se bloquea el mutex nuevamente para realizar operaciones en la matriz dentro de la sección crítica.
            std::unique_lock<std::mutex> lock(mutexAB);
            //Se realiza la operación para hacer ceros por encima y por debajo del pivote en la columna actual.
            factor = AB[posfactor];
        }

        {
            // Bloquear el mutex para garantizar acceso exclusivo a la sección crítica
            std::unique_lock<std::mutex> lock(mutexAB);
            // Realizar operaciones en la matriz dentro de la sección crítica
            for (int j = pospivot; j < posfinfila; ++j) {
                int posactualelim = j + (n + 1) * (id + 1);
                AB[posactualelim] = -1 * factor * AB[j] + AB[posactualelim];
            }
        }
    }
}
```

```c++
/*Parámetros:

AB: Vector que representa la matriz aumentada [A|B].
n: Número de filas de la matriz A.
poscol: Índice de la columna actual en la que se está realizando la eliminación hacia adelante.
Bucle externo:

Itera sobre las filas de la matriz, excluyendo la última fila ya que no es necesario realizar operaciones en la última fila durante la eliminación hacia adelante.
Cálculo de índices:

pospivot: Índice de la posición del pivote en la fila actual.
posfinfila: Índice que marca el final de la fila actual.
Bloqueo del mutex y operaciones en la matriz:

Se bloquea el mutex mutexAB utilizando un std::unique_lock para garantizar acceso exclusivo a la sección crítica.
Se realiza la operación de dividir toda la fila actual por el pivote (piv).
Cálculo de índices y bloqueo del mutex:

posfactor: Índice de la posición del factor en la fila actual.
Se bloquea el mutex nuevamente para obtener el valor del factor.
Bloqueo del mutex y operaciones en la matriz:

Se bloquea el mutex nuevamente para realizar operaciones en la matriz dentro de la sección crítica.
Se realiza la operación para hacer ceros por encima y por debajo del pivote en la columna actual.*/
```

## Funcion EliminacionAtras

La función `eliminacionAtras` realiza la eliminación hacia atrás en el método de Gauss-Jordan. También utiliza un mutex para garantizar el acceso exclusivo a la sección crítica de la matriz AB. Normaliza la última fila y resta múltiplos de esta fila a las filas anteriores.
 Esta fase tiene como objetivo convertir la matriz aumentada [A|B] en una forma escalonada. 
```c++
void eliminacionAtras(std::vector<float>& AB, int n, int poscol) {
    //Itera sobre las filas de la matriz, excluyendo la última fila ya que no es necesario realizar operaciones en la última fila durante la eliminación hacia atrás.
    for (int id = 0; id < n - 1 - poscol; ++id) {
        int pospivot = (n + 2) * (n - 1 - poscol);
        // poscol: Índice de la columna actual en la que se está realizando la eliminación hacia atrás.
        if (poscol == 0) {
            float pivot;

            {
                // Bloquear el mutex para garantizar acceso exclusivo a la sección crítica
                std::unique_lock<std::mutex> lock(mutexAB);
                pivot = AB[pospivot];
                // Realizar operaciones en la matriz dentro de la sección crítica
                AB[pospivot] = AB[pospivot] / pivot;
                AB[pospivot + 1] = AB[pospivot + 1] / pivot;
            }
        }

        float factor; //Valor del elemento en la fila actual por debajo del pivote.
        // indices de las posiciones a las que se les aplicaran las operaciones
        int posactualelim1 = pospivot - (n + 1) * (id + 1); 
        int posactualelim2 = pospivot - (n + 1) * (id + 1) + 1 + poscol;

        {
            // Bloquear el mutex para garantizar acceso exclusivo a la sección crítica
            std::unique_lock<std::mutex> lock(mutexAB);
            factor = AB[pospivot - (n + 1) * (id + 1)];
        }

        {
            // Bloquear el mutex para garantizar acceso exclusivo a la sección crítica
            std::unique_lock<std::mutex> lock(mutexAB);
            // Realizar operaciones en la matriz dentro de la sección crítica
            // estas operaciones están destinadas a hacer ceros por encima y por debajo del pivote en la columna actual
            AB[posactualelim1] = -1 * factor * AB[pospivot] + AB[posactualelim1];
            AB[posactualelim2] = -1 * factor * AB[pospivot + 1 + poscol] + AB[posactualelim2];
        }
    }
}
/*

```


```c++
/*
Parámetros:

AB: Vector que representa la matriz aumentada [A|B].
n: Número de filas de la matriz A.
poscol: Índice de la columna actual en la que se está realizando la eliminación hacia atrás.
Bucle externo:

Itera sobre las filas de la matriz, excluyendo la última fila ya que no es necesario realizar operaciones en la última fila durante la eliminación hacia atrás.
Cálculo de índices:

pospivot: Índice de la posición del pivote en la fila actual.
Condición especial para la primera columna:

Si poscol es igual a cero, se realiza una condición especial para la primera columna:
Se bloquea el mutex mutexAB utilizando un std::unique_lock para garantizar acceso exclusivo a la sección crítica.
Se obtiene el valor del pivote y se realiza la operación de dividir la fila actual por el pivote.
Cálculo de índices y bloqueo del mutex:

factor: Valor del elemento en la fila actual por debajo del pivote.
posactualelim1 y posactualelim2: Índices de las posiciones a las que se les aplicarán las operaciones.
Bloqueo del mutex y operaciones en la matriz:

Se bloquea el mutex utilizando un std::unique_lock para garantizar acceso exclusivo a la sección crítica.
Se realizan operaciones en la matriz dentro de la sección crítica, haciendo ceros por encima y por debajo del pivote en la columna actual.*/
```

## Funcion gaussJordanParallel

la función divide las filas en bloques y ejecuta las fases de eliminación hacia adelante y hacia atrás de manera paralela utilizando std::async. Luego, espera a que todas las tareas paralelas finalicen antes de continuar con la siguiente iteración. Este enfoque paralelo está destinado a mejorar el rendimiento en sistemas multicore al aprovechar la capacidad de ejecución simultánea de varias tareas.
```c++

// Función para realizar el método de Gauss-Jordan de manera paralela
void gaussJordanParallel(std::vector<float>& AB, int n, int numBlocks) {
   // numBlocks: Número de bloques en los que se dividirán las filas para la ejecución paralela.
    // Vector de futuros para almacenar los resultados de las tareas paralelas
    std::vector<std::future<void>> futures;

    // Bucle externo para iterar a través de las columnas excluye la ultima columna porque no es necesario realizar operaciones durante el metodo de gauss jordan ahi
    for (int i = 0; i < n - 1; ++i) {
        // Bucle interno para dividir las filas en bloques y ejecutar tareas en paralelo
        for (int j = 0; j < numBlocks; ++j) {
            // Calcular el rango de filas para el bloque actual
            int start = j * (n / numBlocks);
            int end = (j + 1) * (n / numBlocks);

            // Emplear std::async para lanzar tareas paralelas para eliminación hacia adelante
            futures.emplace_back(std::async(std::launch::async, [&AB, n, i, start, end] {
                for (int id = start; id < end; ++id) {
                    eliminacionAdelante(AB, n, i);
                }
            }));

            // Emplear std::async para lanzar tareas paralelas para eliminación hacia atrás
            futures.emplace_back(std::async(std::launch::async, [&AB, n, i, start, end] {
                for (int id = start; id < end; ++id) {
                    eliminacionAtras(AB, n, i);
                }
            }));
        }

        // Esperar a que todas las tareas paralelas finalicen antes de continuar
        for (auto &future : futures) {
            future.wait();
        }

        // Limpiar el vector de futuros para la próxima iteración
        futures.clear();
    }
}

```
` `
 

## Funcion Secuencial
la función `gaussJordanSecuencial` implementa el método de Gauss-Jordan de manera secuencial, realizando la eliminación hacia adelante y hacia atrás en cada columna de la matriz aumentada [A|B]. A diferencia de la versión paralela, en este caso, las operaciones se realizan de manera secuencial, una fila a la vez, para cada columna.

```c++
// n: numero de filas de la matriz A
void gaussJordanSecuencial(std::vector<float>& AB, int n) {
    for (int i = 0; i < n - 1; ++i) { //bucle externo 
        for (int id = 0; id < n - 1 - i; ++id) { //bucle interno, itera excluyendo las filas ya procesadas anteriormente
            eliminacionAdelante(AB, n, i);
            eliminacionAtras(AB, n, i);
        }
        //Estas funciones realizan las operaciones necesarias para hacer ceros por encima y por debajo del pivote en la columna actual.
    }
}
```

## Funcion Principal 
la función `principal` inicializa matrices, ejecuta el método de Gauss-Jordan de manera secuencial y paralela, mide los tiempos de ejecución y muestra los resultados, incluyendo la aceleración y eficiencia en función del número de hilos utilizados.
```c++

// Función principal
int main() {
    const int n = 100;
    srand(static_cast<unsigned>(time(nullptr))); //se inicializa la semilla del generador de números aleatorios.

    std::vector<float> d_AB, d_AB_secuencial;

    // Inicializar y mostrar la matriz inicial
    inicializarMatriz(d_AB, n + 1, n);
    d_AB_secuencial = d_AB;

    std::cout << "Matriz Inicial (d_AB):" << std::endl;
    imprimirMatriz(d_AB, n + 1, n);

    // Medir el tiempo de ejecución del método de Gauss-Jordan de manera secuencial
    auto startTimeSecuencial = std::chrono::high_resolution_clock::now();
    gaussJordanSecuencial(d_AB_secuencial, n);
    auto endTimeSecuencial = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sequential_time = endTimeSecuencial - startTimeSecuencial; //tiempo transcurrido 

    // Determinar el tamaño del bloque y el número de bloques para la ejecución paralela
    // el mínimo entre el número de hilos disponibles (std::thread::hardware_concurrency()) y un valor fijo de 12.
    int blocksize = std::min(static_cast<int>(std::thread::hardware_concurrency()), 12);
    //Se calcula el número de bloques necesarios para la ejecución paralela basándose en el tamaño de la matriz y el tamaño del bloque.
    int numBlocks = ceil(n / static_cast<float>(blocksize));

    // Mostrar el número de hilos utilizados
    std::cout << "Número de hilos utilizados: " << blocksize << std::endl;

    // Medir el tiempo de ejecución del método de Gauss-Jordan de manera paralela
    auto startTime = std::chrono::high_resolution_clock::now();
    gaussJordanParallel(d_AB, n, numBlocks);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parallel_time = endTime - startTime; //tiempo transcurrido

    // Imprimir resultados
    std::cout << "Tiempo Secuencial: " << sequential_time.count() * 1000 << " ms\n";
    std::cout << "Tiempo Paralelo: " << parallel_time.count() * 1000 << " ms\n";
    std::cout << "Aceleración: " << sequential_time / parallel_time << "\n";
    std::cout << "Eficiencia: " << (sequential_time / parallel_time) / blocksize * 100 << "%\n";

    return 0;
}



```


```c++

```


```c++

```