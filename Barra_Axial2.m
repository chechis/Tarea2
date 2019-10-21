function [k]=Barra_Axial2(E,A,Le)

k=E*A/Le*[1 0 -1 0
          0 0 0 0
          -1 0 1 0
          0 0 0 0];