# Quantum-purity
This repository provides a quantum purity calculation problem for the genetic algorithm provided in https://github.com/vasekp/quantum-ga

The problem is described in detail in https://arxiv.org/pdf/1804.03719.pdf (section "XIV. Quantum Principal Component Analysis").

In short, with this problem the genetic algorithm tries to find a five-qubit quantum circuit that calculates the purity of a one-qubit quantum state. In the procedure, the state of interest |ψ0⟩ is first classically purified to a two-qubit pure state |ψ⟩. As input, the five-qubit circuit takes a joint state of one ancilla qubit in |0⟩ state and two copies of |ψ⟩. After the circuit, the expectation value of the ancilla qubit is measured in the Z-basis, and that expectation value is the purity of the one-qubit state |ψ0⟩. Inputs and outputs are visualized in the following figure:

![Alt Text](inputs_outputs.png)

## The included files

* **QPCA.hpp**     The header file that defines the functions needed for the purity calculation problem.

* **Makefile**     Slightly modified Makefile that recognizes the QPCA.hpp and includes the Armadillo library.

* **quantum.cpp**  Slightly modified main file of the algorithm that recognizes the QPCA.hpp and changes the qubit count.

After succesful installation of the genetic algorithm, this problem works by replacing the Makefile and quantum.cpp files with the ones provided here. The QPCA.hpp file must be placed within the /include/QGA_Problem directory of the genetic algorithm.

##

Adapted from the work in https://github.com/vasekp/quantum-ga by Henrik Romppainen.
