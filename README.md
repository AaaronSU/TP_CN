##########################################################
# README
#
# T. Dufaud
##########################################################

This directory contains the code corresponding to the solution
of Poisson 1D problem by direct method or iterative method.
It is organized in three directories:
src/ 
include/
bin/

"src" contains the source codes, "include" contains the 
header files and "bin" contains the executables. 
The compilation and execution can be done using the Makefile.

Here are the principal targets: 
testenv: bin/tp_testenv
tp2poisson1D_direct: bin/tpPoisson1D_direct
tp2poisson1D_iter: bin/tpPoisson1D_iter

The command,
$ make target
Compile an executable bin/target 

$ make all
compile the executable corresponding to all targets

$ make run_target
Execute ./bin/target

$ make clean
rm *.o bin/*

Here:
$ make run_testenv
$ make run_tpPoisson1D_iter
$ make run_tpPoisson1D_direct


# Rapport sur la résolution de l'équation de la chaleur en 1D stationnaire.


**Matière:** Calcul Numérique  
**Auteur:** Yingqin SU  
**Date:** 8 Janvier 2024

## Introduction

L'équation de chaleur est une équation aux dérivées partielles qui modélise la propagation de la chaleur dans un matériau au fil du temps. Elle joue un rôle crucial dans divers domaines tels que la physique, l'ingénierie, la météorologie et la biologie. Cette équation permet de décrire comment la distribution de la température évolue dans un système donné, en fonction des conditions initiales et des propriétés thermiques du matériau.

Afin d'explorer plus en détail le comportement de l'équation de chaleur, nous allons mettre en œuvre la résolution numérique de l'équation de chaleur unidimensionnelle en utilisant le langage de programmation C. Nous aborderons deux approches distinctes : la méthode directe basée sur la décomposition LU (Lower-Upper) et des méthodes itératives telles que la méthode de Richardson avec un paramètre alpha variable, la méthode de Jacobi ainsi que la méthode de Gauss-Seidel.

Ces approches offrent des perspectives différentes pour résoudre l'équation de chaleur, permettant ainsi une comparaison des performances et une compréhension approfondie de la résolution numérique. À travers cette implémentation, nous explorerons les avantages et les limitations de chaque méthode, contribuant ainsi à une meilleure appréhension de la modélisation thermique dans des contextes variés.

## Equation de la chaleur 1D

L'équation de la chaleur en une dimension, dans le domaine $0 < x < 1$, avec les conditions aux bords $T(0) = T_0$ et $T(1) = T_1$, est formulée comme suit:

$$-k \frac{\partial^2 T}{\partial x^2} = g, \quad \text{pour } x \in ]0,1[$$

où $g$ est le terme source, $k > 0$ est le coefficient de conductivité thermique, $T(x)$ représente la température en fonction de la position $x$, et $T_0$ et $T_1$ sont les températures aux bords du domaine ($T_0 < T_1$).

En notation matricielle, le système discret associé est défini comme:

$$Au = f, \quad A \in \mathbb{R}^{n \times n}, \quad u, f \in \mathbb{R}^n$$

La matrice $A$ est construite à partir des coefficients de différenciation centrée d'ordre 2. Le vecteur $f$ représente les valeurs de température aux nœuds, tandis que le vecteur $u$ contient les termes de source.

La solution analytique de ce problème de la chaleur est donnée par $T(x) = T_0 + x(T_1 - T_0)$.

## Travail préliminaire : Établissement d’un cas de test

En utilisant le développement de Taylor, les expressions suivantes peuvent être déduites :

$$T(x+h) = T(x) + h \frac{dT}{dx} + \frac{h^2}{2} \frac{d^2T}{dx^2}+o(h^2) \quad (1)$$  
$$T(x-h) = T(x) - h \frac{dT}{dx} + \frac{h^2}{2} \frac{d^2T}{dx^2}+o(h^2) \quad (2)$$

En combinant les équations $(1)$ et $(2)$, on obtient :

$$T(x+h) + T(x-h) \approx 2T(x) + h^2 \frac{d^2T}{dx^2}$$

Ainsi, on dérive le résultat suivant :

$$
\begin{align*}
 \frac{d^2T}{dx^2} &\approx \frac{T(x+h) + T(x-h) - 2T(x)}{h^2} \\
                   &\approx \frac{1}{h^2}[[T(x+h) - T(x)] - [T(x)-T(x-h)]]\\
                   &\approx \frac{1}{h^2}[\Delta T(x, h) - \Delta T(x, -h)]\\
                   &\approx \frac{1}{h^2}\Delta^2 T(x, h) 

\end{align*}$$