import os
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import tkinter as tk
import importlib


# Lista de las librerías necesarias
librerias = ["tkinter", "matplotlib"]

# Verificar si las librerías están instaladas
for libreria in librerias:
    try:
        importlib.import_module(libreria)
        print(f"La librería {libreria} está instalada.")
    except ImportError:
        print(f"La librería {libreria} no está instalada. Instalando...")
        try:
            import pip
            pip.main(["install", libreria])
            print(f"La librería {libreria} se instaló correctamente.")
        except Exception as e:
            print(f"No se pudo instalar la librería {libreria}. Error: {e}")


def read_files(carpeta, tipo_grafica):
    datos = {}
    
    for archivo in os.listdir(carpeta):
        if archivo.startswith(tipo_grafica):
            tamano = re.search(r"_(\d+)", archivo).group(1)
            tamano = int(tamano)
            
            with open(os.path.join(carpeta, archivo), "r") as file:
                contenido = file.readlines()
                
            datos_tamano = [float(linea.strip()) for linea in contenido]
            datos[tamano] = datos_tamano
    
    return datos


def build_single_plot(tipo_grafica, N):
    carpeta = "./Clean"
    datos = read_files(carpeta, tipo_grafica)
    
    # Configurar el título y etiquetas de los ejes
    if tipo_grafica == "eficiencia":
        titulo = "Gráfica de Eficiencia"
        etiqueta_eje_y = "Eficiencia"
        grafica_ideal = True
        valor_ideal = 1.0
    elif tipo_grafica == "gflops":
        titulo = "Gráfica de GFLOPS"
        etiqueta_eje_y = "GFLOPS"
        grafica_ideal = False
    elif tipo_grafica == "speedup":
        titulo = "Gráfica de Speedup"
        etiqueta_eje_y = "Speedup"
        grafica_ideal = True
        valor_ideal = None
    elif tipo_grafica == "tiempoCalculo":
        titulo = "Gráfica de Tiempo de Cálculo"
        etiqueta_eje_y = "Tiempo de Cálculo (s)"
        grafica_ideal = False
    elif tipo_grafica == "tiempoComunicacion":
        titulo = "Gráfica de Tiempo de Comunicación"
        etiqueta_eje_y = "Tiempo de Comunicación (s)"
        grafica_ideal = False
    elif tipo_grafica == "tiempoTotal":
        titulo = "Gráfica de Tiempo Total"
        etiqueta_eje_y = "Tiempo Total (s)"
        grafica_ideal = False
    else:
        print("Tipo de gráfica no válido.")
        return
    
    num_procesos = [1, 2, 4, 8, 16, 32, 64]
    
    # Generar las curvas de la gráfica
    colores = cm.get_cmap('viridis', 1)
    formas = ['o', 's', 'D','P', 'H']

    datos_tamano = datos[N]
    curva = []
    for j, proceso in enumerate(num_procesos):
        curva.append(datos_tamano[j])
        
    # Crear la gráfica para el tamaño actual
    plt.plot(num_procesos, curva, '-', label=f"Tamaño {N}", color=colores(0))
    plt.scatter(num_procesos, curva, marker=formas[0], color=colores(0))
    
    # Añadir la gráfica ideal si corresponde
    if grafica_ideal:
        if valor_ideal is not None:
            plt.plot(num_procesos, [valor_ideal] * len(num_procesos), '--', color='red', label="Ideal")
        else:
            plt.plot(num_procesos, num_procesos, '--', color='red', label="Ideal")
    
    # Etiquetas de los ejes y título
    plt.xlabel("Cantidad de Procesos")
    plt.ylabel(etiqueta_eje_y)
    plt.title(titulo)

    plt.legend()
    
    plt.grid(True)
    
    plt.show()

def mostrar_grafica():
    build_single_plot("eficiencia", N.get())

def modificar_N(newValue):
    N.set(newValue)

def calcular_eficiencia():
    carpeta = "./Clean"
    datos = read_files(carpeta, "eficiencia")
    resultado_1_proceso_var.set(datos[N.get()][0])
    resultado_2_procesos_var.set(datos[N.get()][1])
    resultado_4_procesos_var.set(datos[N.get()][2])
    resultado_8_procesos_var.set(datos[N.get()][3])
    resultado_16_procesos_var.set(datos[N.get()][4])
    resultado_32_procesos_var.set(datos[N.get()][5])
    resultado_64_procesos_var.set(datos[N.get()][6])
    
    num_procesos = [1, 2, 4, 8, 16, 32, 64]
    for i in range(len(datos[N.get()])):
        if datos[N.get()][i] < 0.7:
            escalabilidad_var.set(f"Es escalable hasta {num_procesos[i-1]} procesos")
            break

# Crear la ventana del HMI
ventana = tk.Tk()
ventana.title("Simulación N-Body")
ventana.geometry("500x400")

# Crear un Frame para el N
N_frame = tk.LabelFrame(ventana, text="Selecciona el N", width=120, height=250)
N_frame.pack(side=tk.LEFT, padx=40, pady=10)
N_frame.pack_propagate(False)

initial_N = 4096
N = tk.IntVar()
N.set(initial_N)
valor_actual_N = tk.Label(N_frame, textvariable=N)
valor_actual_N.pack(pady=10)

# Botón para N = 4096
boton_4096_N = tk.Button(N_frame, text="4096", command=lambda: modificar_N(4096))
boton_4096_N.pack(pady=5)

# Botón para N = 8192
boton_8192_N = tk.Button(N_frame, text="8192", command=lambda: modificar_N(8192))
boton_8192_N.pack(pady=5)

# Botón para N = 16384
boton_16384_N = tk.Button(N_frame, text="16384", command=lambda: modificar_N(16384))
boton_16384_N.pack(pady=5)

# Botón para N = 32768
boton_32768_N = tk.Button(N_frame, text="32768", command=lambda: modificar_N(32768))
boton_32768_N.pack(pady=5)

# Botón para N = 65536
boton_65536_N = tk.Button(N_frame, text="65536", command=lambda: modificar_N(65536))
boton_65536_N.pack(pady=10)


# Botón para calcular los resultados
boton_calcular = tk.Button(ventana, text="Calcular", command=calcular_eficiencia)
boton_calcular.place(relx=0.33, rely=0.9) 

# Botón para calcular los resultados
boton_graficar = tk.Button(ventana, text="Gráfica", command=mostrar_grafica)
boton_graficar.place(relx=0.66, rely=0.9) 

# Crear un Frame para los resultados
resultados_frame = tk.LabelFrame(ventana, text="Eficiencia por proceso de menor a mayor")
resultados_frame.pack(side=tk.LEFT, padx=10, pady=10)

# Variables para los resultados
resultado_1_proceso_var = tk.IntVar()
resultado_2_procesos_var = tk.IntVar()
resultado_4_procesos_var = tk.IntVar()
resultado_8_procesos_var = tk.IntVar()
resultado_16_procesos_var = tk.IntVar()
resultado_32_procesos_var = tk.IntVar()
resultado_64_procesos_var = tk.IntVar()
escalabilidad_var = tk.StringVar()

resultado_1_proceso = tk.Label(resultados_frame, textvariable=resultado_1_proceso_var)
resultado_1_proceso.pack()

resultado_2_procesos = tk.Label(resultados_frame, textvariable=resultado_2_procesos_var)
resultado_2_procesos.pack()

resultado_4_procesos = tk.Label(resultados_frame, textvariable=resultado_4_procesos_var)
resultado_4_procesos.pack()

resultado_8_procesos = tk.Label(resultados_frame, textvariable=resultado_8_procesos_var)
resultado_8_procesos.pack()

resultado_16_procesos = tk.Label(resultados_frame, textvariable=resultado_16_procesos_var)
resultado_16_procesos.pack()

resultado_32_procesos = tk.Label(resultados_frame, textvariable=resultado_32_procesos_var)
resultado_32_procesos.pack()

resultado_64_procesos = tk.Label(resultados_frame, textvariable=resultado_64_procesos_var)
resultado_64_procesos.pack()

escalabilidad = tk.Label(resultados_frame, textvariable=escalabilidad_var)
escalabilidad.pack()

# Iniciar la ventana
ventana.mainloop()