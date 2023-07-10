# Evolucion_Galactica_Paralela

[Link Overleaf](https://www.overleaf.com/9284392586hypdznkwnfhp)

## Secuencia de pasos para ejecutar en Khipu
1. Git clone repository
2. cd Evolucion_Galactica_Paralela
3. cd N-Body
4. rm cpu-4th
5. module load gcc/5.5.0
6. module load mpich/4.0
7. (opcional editar archivo para cambiar el N) nano gen-plum.c
8. gcc -o new_n gen-plum.c -lm
9. ./new_n
10. make cpu-4th
11. mpirun -np 2 ./cpu-4th