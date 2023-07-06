import os
import re

carpeta_lectura = "./Raw"
carpeta_escritura = "./Clean"

# Diccionarios para almacenar los datos por tamaño n del problema
datos_speedup = {}
datos_eficiencia = {}

# Leer los archivos de la carpeta en orden de tamaño de procesos
archivos = sorted(os.listdir(carpeta_lectura), key=lambda x: (int(re.findall(r"(\d+)", x)[0]), int(re.findall(r"(\d+)", x)[1])))
for archivo in archivos:
    if archivo.endswith(".txt"):
        # Obtener el tamaño y el número de procesos del archivo
        tamano, num_procesos = re.findall(r"(\d+)", archivo)
        tamano = int(tamano)
        num_procesos = int(num_procesos)
        #print(tamano, num_procesos)
        
        # Leer el contenido del archivo
        with open(os.path.join(carpeta_lectura, archivo), "r") as file:
            contenido = file.read()

        # Extraer los datos necesarios
        real_speed = float(re.search(r"Real Speed = ([\d.]+) GFlops", contenido).group(1))
        tiempo_tot = float(re.search(r"tot  :  ([\d.E+-]+)", contenido).group(1))
        tiempo_comm = float(re.search(r"comm :  ([\d.E+-]+)", contenido).group(1))
        tiempo_force = float(re.search(r"force:  ([\d.E+-]+)", contenido).group(1))

        if num_procesos == 1:
            tiempo_secuencial = tiempo_tot 
        # Calcular el speedup y la eficiencia
        speedup = tiempo_secuencial / tiempo_tot
        eficiencia = speedup / num_procesos
        
        # Guardar speedup y eficiencia en diccionarios
        if tamano not in datos_speedup:
            datos_speedup[tamano] = []
        if tamano not in datos_eficiencia:
            datos_eficiencia[tamano] = []
        
        datos_speedup[tamano].append(speedup)
        datos_eficiencia[tamano].append(eficiencia)
        
        # Guardar los datos adicionales en archivos
        archivo_gflops = os.path.join(carpeta_escritura, f"gflops_{tamano}.txt")
        with open(archivo_gflops, "a") as file:
            file.write(f"{real_speed:.8f}\n")
        
        archivo_tiempo_total = os.path.join(carpeta_escritura, f"tiempoTotal_{tamano}.txt")
        with open(archivo_tiempo_total, "a") as file:
            file.write(f"{tiempo_tot:.8f}\n")
        
        archivo_tiempo_comm = os.path.join(carpeta_escritura, f"tiempoComunicacion_{tamano}.txt")
        with open(archivo_tiempo_comm, "a") as file:
            file.write(f"{tiempo_comm:.8f}\n")
        
        archivo_tiempo_calc = os.path.join(carpeta_escritura, f"tiempoCalculo_{tamano}.txt")
        with open(archivo_tiempo_calc, "a") as file:
            file.write(f"{tiempo_force:.8f}\n")
            
for tamano, speedup_lista in datos_speedup.items():
    archivo_speedup = os.path.join(carpeta_escritura, f"speedup_{tamano}.txt")
    with open(archivo_speedup, "w") as file:
        for speedup in speedup_lista:
            file.write(f"{speedup:.8f}\n")

for tamano, eficiencia_lista in datos_eficiencia.items():
    archivo_eficiencia = os.path.join(carpeta_escritura, f"eficiencia_{tamano}.txt")
    with open(archivo_eficiencia, "w") as file:
        for eficiencia in eficiencia_lista:
            file.write(f"{eficiencia:.8f}\n")
