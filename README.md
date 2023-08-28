# Simulación de Sistemas de Comunicaciones Inalámbricas MIMO
En este repositorio se encuentran rutinas de código necesarias para emular el comportamiento de sistemas de comunicaciones inalámbricas MIMO en el lenguaje de programación Python. Particularmente se obtiene la tasa de error de bit promedio (ABEP) vs relación señal a ruido (SNR) en decibeles.
Se compara el rendimiento de manera analítica y simulada en cada uno de los métodos, en los cuales se tienen diversos sistemas al variar el número de antenas en el transmisor (Tx) y/o receptor (Rx).
Las librerías de Python necesarias para ejecutar las rutinas son:
* Numpy
* Matplotlib
* Scipy

A continuación se describen las rutinas:

* 4QAM_SISO.py
  * Sistema de 1 antena en el transmisor y 1 en el receptor. Se emplea una modulación 4-QAM

*  16QAM_SISO.py
   * Sistema de 1 antena en el transmisor y 1 en el receptor. Se emplea una modulación 16-QAM

* 4QAM_SIMO_1x2.py
  * Sistema de 1 antena en el transmisor y 2 antenas en el receptor. Se emplea una modulación 4-QAM.

* 4QAM_MISO_2x1.py
  * Sistema de 2 antenas en el transmisor y 1 antena en el receptor. Se emplea una modulación 4-QAM.
 
* 4QAM_MIMO_2x2.py
  * Sistema de 2 antenas en el transmisor y 2 en el receptor. Se emplea una modulación 4-QAM.

