import os
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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

def build_plots(carpeta, tipo_grafica):
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
    
    tamanos = sorted(datos.keys())
    num_procesos = [1, 2, 4, 8, 16, 32, 64]
    
    # Generar las curvas de la gráfica
    colores = cm.get_cmap('viridis', len(tamanos))
    formas = ['o', 's', 'D','P', 'H']
    for i, tamano in enumerate(tamanos):
        datos_tamano = datos[tamano]
        curva = []
        for j, proceso in enumerate(num_procesos):
            curva.append(datos_tamano[j])
        
        # Crear la gráfica para el tamaño actual
        plt.plot(num_procesos, curva, '-'+formas[i%len(formas)], label=f"Tamaño {tamano}", color=colores(i))
        plt.scatter(num_procesos, curva, marker=formas[i%len(formas)], color=colores(i))
    
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

# Ejemplo de uso
carpeta_datos = "./Clean"
tipo_grafica = "eficiencia"

build_plots(carpeta_datos, tipo_grafica)
