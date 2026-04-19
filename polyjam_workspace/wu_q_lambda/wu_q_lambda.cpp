#include <polyjam/polyjam.hpp>

namespace
{

using namespace polyjam;
using namespace polyjam::core;
using std::list;
using std::string;
using std::stringstream;
using std::vector;

const int kQuaternionSize = 4;
const int kCubicMonomialCount = 64;
const int kWRows = 4;
const int kWCols = 64;
const int kQRows = 4;
const int kQCols = 4;
const int kCoeffCount = kWRows * kWCols + kQRows * kQCols;

inline int cubicIndex(int a, int b, int c)
{
  return (a * kQuaternionSize + b) * kQuaternionSize + c;
}

Poly nonzeroRandZ(size_t nu)
{
  for(int attempt = 0; attempt < 1024; ++attempt)
  {
    Poly value = Poly::randZ(nu);
    if(!value.isZero())
      return value;
  }
  return Poly::oneZ(nu);
}

PolyMatrix cubicQuaternionMonomials(PolyMatrix & q)
{
  const size_t nu = q[0].leadingTerm().monomial().dimensions();
  PolyMatrix u(Poly::zeroSZ(nu),kCubicMonomialCount,1);

  for(int a = 0; a < kQuaternionSize; ++a)
    for(int b = 0; b < kQuaternionSize; ++b)
      for(int c = 0; c < kQuaternionSize; ++c)
        u[cubicIndex(a,b,c)] = q[a] * q[b] * q[c];

  return u;
}

PolyMatrix cubicQuaternionMonomialsZ(PolyMatrix & q)
{
  const size_t nu = q[0].leadingTerm().monomial().dimensions();
  PolyMatrix u(Poly::zeroZ(nu),kCubicMonomialCount,1);

  for(int a = 0; a < kQuaternionSize; ++a)
    for(int b = 0; b < kQuaternionSize; ++b)
      for(int c = 0; c < kQuaternionSize; ++c)
        u[cubicIndex(a,b,c)] = q[a] * q[b] * q[c];

  return u;
}

PolyMatrix randomUnitQuaternionZ(size_t nu)
{
  // Rational stereographic parametrization of S^3:
  //   q = [(1-r^2), 2x, 2y, 2z] / (1+r^2).
  // It avoids finite-field square roots and exactly satisfies q^T q = 1.
  for(int attempt = 0; attempt < 2048; ++attempt)
  {
    Poly x = Poly::randZ(nu);
    Poly y = Poly::randZ(nu);
    Poly z = Poly::randZ(nu);

    Poly r2 = x*x + y*y + z*z;
    Poly denom = Poly::oneZ(nu) + r2;
    if(denom.isZero())
      continue;

    Poly invDenom = Poly::oneZ(nu).leadingTerm() / denom.leadingTerm();

    PolyMatrix q(Poly::zeroZ(nu),kQuaternionSize,1);
    q[0] = (Poly::oneZ(nu) - r2) * invDenom;
    q[1] = Poly::constZ(2,nu) * x * invDenom;
    q[2] = Poly::constZ(2,nu) * y * invDenom;
    q[3] = Poly::constZ(2,nu) * z * invDenom;

    // q0 is used for the synthetic Q construction below, and qz is the
    // default action variable used by Polyjam.  Avoid degenerate samples.
    if(!q[0].isZero() && !q[3].isZero())
      return q;
  }

  PolyMatrix q(Poly::zeroZ(nu),kQuaternionSize,1);
  q[0] = Poly::oneZ(nu);
  return q;
}

bool allCoefficientsNonzero(const vector<Poly> & coeffs)
{
  if(coeffs.size() != static_cast<size_t>(kCoeffCount))
    return false;

  for(size_t i = 0; i < coeffs.size(); ++i)
    if(coeffs[i].isZero())
      return false;

  return true;
}

bool buildConsistentCoefficientSample(size_t nu, vector<Poly> & coeffs)
{
  coeffs.clear();
  coeffs.resize(kCoeffCount,Poly::zeroZ(nu));

  const int wOffset = 0;
  const int qOffset = kWRows * kWCols;

  for(int attempt = 0; attempt < 2048; ++attempt)
  {
    PolyMatrix qgt = randomUnitQuaternionZ(nu);
    if(qgt[0].isZero())
      continue;

    Poly lambda = nonzeroRandZ(nu);
    PolyMatrix ugt = cubicQuaternionMonomialsZ(qgt);

    vector<Poly> candidate(kCoeffCount,Poly::zeroZ(nu));

    // Fill W with non-zero finite-field values.
    for(int r = 0; r < kWRows; ++r)
      for(int c = 0; c < kWCols; ++c)
        candidate[wOffset + r*kWCols + c] = nonzeroRandZ(nu);

    // Fill columns 1..3 of Q randomly.  Column 0 is then set so that the
    // selected ground-truth (lambda,q) exactly satisfies
    // W*u(q) = (Q-lambda*I)*q.
    for(int r = 0; r < kQRows; ++r)
    {
      Poly target = Poly::zeroZ(nu);
      for(int c = 0; c < kWCols; ++c)
        target += candidate[wOffset + r*kWCols + c] * ugt[c];

      target += lambda * qgt[r];

      Poly known = Poly::zeroZ(nu);
      for(int c = 1; c < kQCols; ++c)
      {
        Poly qrc = nonzeroRandZ(nu);
        candidate[qOffset + r*kQCols + c] = qrc;
        known += qrc * qgt[c];
      }

      Poly invQ0 = Poly::oneZ(nu).leadingTerm() / qgt[0].leadingTerm();
      candidate[qOffset + r*kQCols + 0] = (target - known) * invQ0;
    }

    if(!allCoefficientsNonzero(candidate))
      continue;

    coeffs.swap(candidate);
    return true;
  }

  coeffs.clear();
  return false;
}

Poly liftCoefficient(const Poly & zpValue, size_t index, size_t nu)
{
  stringstream name;
  name << "coeffs[" << index << "]";
  return Poly(Term(
      Coefficient(name.str()),
      zpValue.leadingTerm().coefficient().clone(),
      Monomial(nu)));
}

} // namespace

int main(int argc, char ** argv)
{
  (void)argc;
  (void)argv;

  initGenerator();

  // Unknown order in generated solutions:
  //   x1 = lambda, x2..x5 = [q_w,q_x,q_y,q_z].
  const size_t nu = 5;

  std::cout << "Generating Groebner-basis solver for W*u(q) = "
            << "(Q - lambda*I_4)*q with u=q kron q kron q." << std::endl;

  vector<Poly> coeffSample;
  if(!buildConsistentCoefficientSample(nu,coeffSample))
  {
    std::cerr << "Unable to build a non-degenerate finite-field sample for "
              << "the W*u=(Q-lambda*I)q system." << std::endl;
    return 1;
  }

  std::cout << "Using " << coeffSample.size()
            << " scalar coefficients: 256 for W and 16 for Q." << std::endl;

  Poly lambda = Poly::uSZ(1,nu);
  PolyMatrix q(Poly::zeroSZ(nu),kQuaternionSize,1);
  q[0] = Poly::uSZ(2,nu);
  q[1] = Poly::uSZ(3,nu);
  q[2] = Poly::uSZ(4,nu);
  q[3] = Poly::uSZ(5,nu);

  PolyMatrix u = cubicQuaternionMonomials(q);

  list<Poly*> eqs;

  const int wOffset = 0;
  const int qOffset = kWRows * kWCols;

  for(int r = 0; r < kWRows; ++r)
  {
    Poly residual = Poly::zeroSZ(nu);

    for(int c = 0; c < kWCols; ++c)
    {
      const int coeffIndex = wOffset + r*kWCols + c;
      residual += liftCoefficient(coeffSample[coeffIndex],coeffIndex,nu) * u[c];
    }

    for(int c = 0; c < kQCols; ++c)
    {
      const int coeffIndex = qOffset + r*kQCols + c;
      residual -= liftCoefficient(coeffSample[coeffIndex],coeffIndex,nu) * q[c];
    }

    residual += lambda * q[r];
    eqs.push_back(new Poly(residual));
  }

  Poly unitQuaternion = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] -
                        Poly::oneSZ(nu);
  eqs.push_back(new Poly(unitQuaternion));

  std::cout << "Constructed " << eqs.size() << " polynomial equations in "
            << nu << " unknowns: lambda and a unit quaternion." << std::endl;

  const string parameters("const std::vector<double> & coeffs");
  execGenerator(eqs,string("wu_q_lambda_core"),parameters,false);

  for(list<Poly*>::iterator it = eqs.begin(); it != eqs.end(); ++it)
    delete *it;

  return 0;
}
