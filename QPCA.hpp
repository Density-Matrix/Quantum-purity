// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

#include "QIClib"

namespace {

using QGA::Backend::State;

using Gene = QGA::Gene<QGA::Gates::CNOT::WithControls<QGA::Controls::ANY>, QGA::Gates::XYZ::WithControls<QGA::Controls::ANY>, QGA::Gates::SWAP, QGA::Gates::Fixed::WithControls<QGA::Controls::ANY>>;
class Candidate : public QGA::CandidateBase<Candidate, Gene, double, unsigned>
{

  using Base = QGA::CandidateBase<Candidate, Gene, double, unsigned>;

public:

  using Base::Base;

  Base::Fitness fitness() const {
    if(genotype().size() > 1000)
      return {};

    double total_error = 0; //Total error of the circuit
    unsigned dim = 1 << Config::nBit; //*** Works correctly only with nBit == 5 ***
    unsigned rho_size = (Config::nBit-1)/2; //The size of one copy of the random state (2x2)
    arma::cx_mat identity = qic::spm.S(0); //Identity matrix
    arma::cx_mat Z1 = qic::spm.S(3); //Pauli sigma_z matrix
    for (unsigned i = 1; i < Config::nBit; i++) { //The (32x32) matrix that measures observable Z on the first qubit only
      Z1 = arma::kron(Z1,identity);
    }
    arma::cx_vec ancilla = qic::spm.basis2(0, 0); //The ancilla qubit |0>
    unsigned checks = 100; //Number of random states to test the circuit with
    arma::vec eigvals(rho_size); //For storing eigenvalues
    arma::cx_mat eigvecs(rho_size,rho_size); //For storing eigenvectors
    arma::cx_vec psi(rho_size*rho_size); //For the purified state
    State in{}; //The state to use as input
    arma::cx_vec arma_out(dim); //The output state in Armadillo format

    //This loop tests the circuit with random states [checks] times
    for(unsigned i = 0; i < checks; i++) {
      // generate a random density matrix rho
      arma::cx_mat rho = qic::randRho(rho_size);
      // calculate a purified state psi (2 qubits) from rho (1 qubit)
      arma::eig_sym(eigvals, eigvecs, rho);
      psi.zeros();
      for (unsigned i = 0; i < rho_size; i++) {
        psi += sqrt(eigvals(i))*(arma::kron(eigvecs.col(i),eigvecs.col(i)));
      }
      // calculate the correct purity from rho
      double correct = std::abs(arma::trace(rho*rho));
      // run a |0> ancilla + 2 copies of psi through the circuit
      arma::cx_vec arma_in = arma::kron(ancilla,arma::kron(psi,psi));
      for (unsigned i = 0; i < dim; i++) { //Convert arma vector to QGA vector
        in[i] = arma_in(i);
      }
      State out = sim(in);
      for (unsigned i = 0; i < dim; i++) { //Convert QGA vector to arma vector
        arma_out(i) = out[i];
      }
      // calculate the purity as the expectation value <psi|Z1|psi>
      double purity = std::abs(arma::as_scalar(arma::trans(arma_out)*Z1*arma_out));  
      // calculate error of purity and add it to total error
      double error = std::abs((purity - correct) / correct);
      total_error += error;
    }
    double average_error = total_error / checks;
    return {
      trimError(average_error),
      genotype().size()
    };
  }

  std::ostream& print_full(std::ostream& os) const {
    os << '\n';
    os << "print_full unimplemented";
    os << '\n';
    return os;
  }

private:

  State sim(const State& psi) const {
    State ret{psi};
    for(const auto& g : genotype())
      ret = g->applyTo(ret);
    return ret;
  }

}; // class Candidate

} // anonymous namespace

#endif // !defined QGA_PROBLEM_HPP
