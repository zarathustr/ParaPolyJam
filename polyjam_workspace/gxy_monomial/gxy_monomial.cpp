#include <polyjam/polyjam.hpp>

namespace
{

using namespace polyjam;
using namespace polyjam::core;
using std::list;
using std::string;
using std::stringstream;
using std::vector;

const int kUnknownCount = 6;
const int kM1Rows = 3;
const int kM1Cols = 40;
const int kM2Rows = 3;
const int kM2Cols = 56;
const int kCoeffCount = kM1Rows * kM1Cols + kM2Rows * kM2Cols;

inline int m1CoeffIndex(int row, int col)
{
  return row * kM1Cols + col;
}

inline int m2CoeffIndex(int row, int col)
{
  return kM1Rows * kM1Cols + row * kM2Cols + col;
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

PolyMatrix randomNonzeroUnknownsZ(size_t nu)
{
  PolyMatrix g(Poly::zeroZ(nu),kUnknownCount,1);
  for(int i = 0; i < kUnknownCount; ++i)
    g[i] = nonzeroRandZ(nu);
  return g;
}

vector<Poly> mons1(PolyMatrix & g)
{
  vector<Poly> m;
  m.reserve(kM1Cols);

  Poly & gX1 = g[0];
  Poly & gX2 = g[1];
  Poly & gX3 = g[2];
  Poly & gY1 = g[3];
  Poly & gY2 = g[4];
  Poly & gY3 = g[5];

  m.push_back(gX1*gY1*gY1);
  m.push_back(gX1*gY1*gY2);
  m.push_back(gX1*gY1*gY3);
  m.push_back(gX1*gY1);
  m.push_back(gX1*gY2*gY2);
  m.push_back(gX1*gY2*gY3);
  m.push_back(gX1*gY2);
  m.push_back(gX1*gY3*gY3);
  m.push_back(gX1*gY3);
  m.push_back(gX1);
  m.push_back(gX2*gY1*gY1);
  m.push_back(gX2*gY1*gY2);
  m.push_back(gX2*gY1*gY3);
  m.push_back(gX2*gY1);
  m.push_back(gX2*gY2*gY2);
  m.push_back(gX2*gY2*gY3);
  m.push_back(gX2*gY2);
  m.push_back(gX2*gY3*gY3);
  m.push_back(gX2*gY3);
  m.push_back(gX2);
  m.push_back(gX3*gY1*gY1);
  m.push_back(gX3*gY1*gY2);
  m.push_back(gX3*gY1*gY3);
  m.push_back(gX3*gY1);
  m.push_back(gX3*gY2*gY2);
  m.push_back(gX3*gY2*gY3);
  m.push_back(gX3*gY2);
  m.push_back(gX3*gY3*gY3);
  m.push_back(gX3*gY3);
  m.push_back(gX3);
  m.push_back(gY1*gY1);
  m.push_back(gY1*gY2);
  m.push_back(gY1*gY3);
  m.push_back(gY1);
  m.push_back(gY2*gY2);
  m.push_back(gY2*gY3);
  m.push_back(gY2);
  m.push_back(gY3*gY3);
  m.push_back(gY3);
  m.push_back(g[0].one());

  return m;
}

vector<Poly> mons1Z(PolyMatrix & g)
{
  vector<Poly> m = mons1(g);
  return m;
}

vector<Poly> mons4(PolyMatrix & g)
{
  vector<Poly> m;
  m.reserve(kM2Cols);

  Poly & gX1 = g[0];
  Poly & gX2 = g[1];
  Poly & gX3 = g[2];
  Poly & gY1 = g[3];
  Poly & gY2 = g[4];
  Poly & gY3 = g[5];

  m.push_back(gX1*gX1*gY1);
  m.push_back(gX1*gX1*gY2);
  m.push_back(gX1*gX1*gY3);
  m.push_back(gX1*gX1);
  m.push_back(gX1*gX2*gY1);
  m.push_back(gX1*gX2*gY2);
  m.push_back(gX1*gX2*gY3);
  m.push_back(gX1*gX2);
  m.push_back(gX1*gX3*gY1);
  m.push_back(gX1*gX3*gY2);
  m.push_back(gX1*gX3*gY3);
  m.push_back(gX1*gX3);
  m.push_back(gX1*gY1);
  m.push_back(gX1*gY2);
  m.push_back(gX1*gY3);
  m.push_back(gX1);
  m.push_back(gX2*gX2*gY1);
  m.push_back(gX2*gX2*gY2);
  m.push_back(gX2*gX2*gY3);
  m.push_back(gX2*gX2);
  m.push_back(gX2*gX3*gY1);
  m.push_back(gX2*gX3*gY2);
  m.push_back(gX2*gX3*gY3);
  m.push_back(gX2*gX3);
  m.push_back(gX2*gY1);
  m.push_back(gX2*gY2);
  m.push_back(gX2*gY3);
  m.push_back(gX2);
  m.push_back(gX3*gX3*gY1);
  m.push_back(gX3*gX3*gY2);
  m.push_back(gX3*gX3*gY3);
  m.push_back(gX3*gX3);
  m.push_back(gX3*gY1);
  m.push_back(gX3*gY2);
  m.push_back(gX3*gY3);
  m.push_back(gX3);
  m.push_back(gY1*gY1*gY1);
  m.push_back(gY1*gY1*gY2);
  m.push_back(gY1*gY1*gY3);
  m.push_back(gY1*gY1);
  m.push_back(gY1*gY2*gY2);
  m.push_back(gY1*gY2*gY3);
  m.push_back(gY1*gY2);
  m.push_back(gY1*gY3*gY3);
  m.push_back(gY1*gY3);
  m.push_back(gY1);
  m.push_back(gY2*gY2*gY2);
  m.push_back(gY2*gY2*gY3);
  m.push_back(gY2*gY2);
  m.push_back(gY2*gY3*gY3);
  m.push_back(gY2*gY3);
  m.push_back(gY2);
  m.push_back(gY3*gY3*gY3);
  m.push_back(gY3*gY3);
  m.push_back(gY3);
  m.push_back(g[0].one());

  return m;
}

vector<Poly> mons4Z(PolyMatrix & g)
{
  vector<Poly> m = mons4(g);
  return m;
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

  for(int attempt = 0; attempt < 4096; ++attempt)
  {
    PolyMatrix ggt = randomNonzeroUnknownsZ(nu);
    vector<Poly> m1gt = mons1Z(ggt);
    vector<Poly> m2gt = mons4Z(ggt);

    vector<Poly> candidate(kCoeffCount,Poly::zeroZ(nu));

    for(int r = 0; r < kM1Rows; ++r)
    {
      Poly sum = Poly::zeroZ(nu);
      for(int c = 0; c < kM1Cols - 1; ++c)
      {
        candidate[m1CoeffIndex(r,c)] = nonzeroRandZ(nu);
        sum += candidate[m1CoeffIndex(r,c)] * m1gt[c];
      }
      candidate[m1CoeffIndex(r,kM1Cols-1)] = Poly::zeroZ(nu) - sum;
    }

    for(int r = 0; r < kM2Rows; ++r)
    {
      Poly sum = Poly::zeroZ(nu);
      for(int c = 0; c < kM2Cols - 1; ++c)
      {
        candidate[m2CoeffIndex(r,c)] = nonzeroRandZ(nu);
        sum += candidate[m2CoeffIndex(r,c)] * m2gt[c];
      }
      candidate[m2CoeffIndex(r,kM2Cols-1)] = Poly::zeroZ(nu) - sum;
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
  //   x1..x6 = [gX1,gX2,gX3,gY1,gY2,gY3].
  const size_t nu = kUnknownCount;

  std::cout << "Generating Groebner-basis solver for the MATLAB M1/M2 "
            << "monomial system in [gX1,gX2,gX3,gY1,gY2,gY3]." << std::endl;

  vector<Poly> coeffSample;
  if(!buildConsistentCoefficientSample(nu,coeffSample))
  {
    std::cerr << "Unable to build a non-degenerate finite-field sample for "
              << "the M1/M2 monomial system." << std::endl;
    return 1;
  }

  std::cout << "Using " << coeffSample.size()
            << " scalar coefficients: 120 for M1 and 168 for M2." << std::endl;

  PolyMatrix g(Poly::zeroSZ(nu),kUnknownCount,1);
  g[0] = Poly::uSZ(1,nu);
  g[1] = Poly::uSZ(2,nu);
  g[2] = Poly::uSZ(3,nu);
  g[3] = Poly::uSZ(4,nu);
  g[4] = Poly::uSZ(5,nu);
  g[5] = Poly::uSZ(6,nu);

  vector<Poly> m1 = mons1(g);
  vector<Poly> m2 = mons4(g);

  list<Poly*> eqs;

  for(int r = 0; r < kM1Rows; ++r)
  {
    Poly eq = Poly::zeroSZ(nu);
    for(int c = 0; c < kM1Cols; ++c)
    {
      const int idx = m1CoeffIndex(r,c);
      eq += liftCoefficient(coeffSample[idx],idx,nu) * m1[c];
    }
    eqs.push_back(new Poly(eq));
  }

  for(int r = 0; r < kM2Rows; ++r)
  {
    Poly eq = Poly::zeroSZ(nu);
    for(int c = 0; c < kM2Cols; ++c)
    {
      const int idx = m2CoeffIndex(r,c);
      eq += liftCoefficient(coeffSample[idx],idx,nu) * m2[c];
    }
    eqs.push_back(new Poly(eq));
  }

  std::cout << "Constructed " << eqs.size() << " polynomial equations in "
            << nu << " unknowns." << std::endl;

  const string parameters("const std::vector<double> & coeffs");
  execGenerator(eqs,string("gxy_monomial_core"),parameters,false);

  for(list<Poly*>::iterator it = eqs.begin(); it != eqs.end(); ++it)
    delete *it;

  return 0;
}
