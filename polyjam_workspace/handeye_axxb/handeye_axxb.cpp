#include <polyjam/polyjam.hpp>

#ifndef POLYJAM_HAND_EYE_PAIRS
#define POLYJAM_HAND_EYE_PAIRS 2
#endif

#if POLYJAM_HAND_EYE_PAIRS < 2
#error "POLYJAM_HAND_EYE_PAIRS must be at least 2 for a zero-dimensional AX=XB hand-eye problem."
#endif

namespace
{

using namespace polyjam;
using namespace polyjam::core;
using std::string;
using std::stringstream;
using std::vector;

inline int ridx(int r, int c) { return 3*r + c; }

PolyMatrix quaternionRotation(
    const Poly & qw,
    const Poly & qx,
    const Poly & qy,
    const Poly & qz)
{
  const size_t nu = qw.leadingTerm().monomial().dimensions();
  PolyMatrix R(qw.zero(),3,3);

  const Poly one = qw.one();
  const Poly two = Poly::constSZ(2,nu);

  R(0,0) = one - two * (qy*qy + qz*qz);
  R(0,1) = two * (qx*qy - qw*qz);
  R(0,2) = two * (qx*qz + qw*qy);

  R(1,0) = two * (qx*qy + qw*qz);
  R(1,1) = one - two * (qx*qx + qz*qz);
  R(1,2) = two * (qy*qz - qw*qx);

  R(2,0) = two * (qx*qz - qw*qy);
  R(2,1) = two * (qy*qz + qw*qx);
  R(2,2) = one - two * (qx*qx + qy*qy);

  return R;
}

PolyMatrix quaternionRotationZ(
    const Poly & qw,
    const Poly & qx,
    const Poly & qy,
    const Poly & qz)
{
  const size_t nu = qw.leadingTerm().monomial().dimensions();
  PolyMatrix R(qw.zero(),3,3);

  const Poly one = qw.one();
  const Poly two = Poly::constZ(2,nu);

  R(0,0) = one - two * (qy*qy + qz*qz);
  R(0,1) = two * (qx*qy - qw*qz);
  R(0,2) = two * (qx*qz + qw*qy);

  R(1,0) = two * (qx*qy + qw*qz);
  R(1,1) = one - two * (qx*qx + qz*qz);
  R(1,2) = two * (qy*qz - qw*qx);

  R(2,0) = two * (qx*qz - qw*qy);
  R(2,1) = two * (qy*qz + qw*qx);
  R(2,2) = one - two * (qx*qx + qy*qy);

  return R;
}

PolyMatrix randomUnitQuaternionZ(size_t nu)
{
  // Rational stereographic parametrization of S^3:
  //   q = [(1-r^2), 2u, 2v, 2w] / (1+r^2).
  // This avoids finite-field square roots while satisfying q^T q = 1 exactly.
  for( int attempt = 0; attempt < 1024; ++attempt )
  {
    Poly u = Poly::randZ(nu);
    Poly v = Poly::randZ(nu);
    Poly w = Poly::randZ(nu);

    // qz is used as the action variable. Avoid a degenerate random instance.
    if( w.isZero() )
      continue;

    Poly r2 = u*u + v*v + w*w;
    Poly denom = Poly::oneZ(nu) + r2;
    if( denom.isZero() )
      continue;

    Poly denomInv = Poly::oneZ(nu).leadingTerm() / denom.leadingTerm();

    PolyMatrix q(Poly::zeroZ(nu),4,1);
    q[0] = (Poly::oneZ(nu) - r2) * denomInv;
    q[1] = Poly::constZ(2,nu) * u * denomInv;
    q[2] = Poly::constZ(2,nu) * v * denomInv;
    q[3] = Poly::constZ(2,nu) * w * denomInv;
    return q;
  }

  // Extremely unlikely fallback for the default prime field.
  PolyMatrix q(Poly::zeroZ(nu),4,1);
  q[0] = Poly::oneZ(nu);
  return q;
}

PolyMatrix randomTranslationZ(size_t nu)
{
  PolyMatrix t(Poly::zeroZ(nu),3,1);
  for( int i = 0; i < 3; ++i )
    t[i] = Poly::randZ(nu);
  return t;
}

PolyMatrix makeTransform(
    PolyMatrix & R,
    PolyMatrix & t,
    size_t nu)
{
  PolyMatrix T(R[0].zero(),4,4);
  for( int r = 0; r < 3; ++r )
  {
    for( int c = 0; c < 3; ++c )
      T(r,c) = R(r,c);
    T(r,3) = t[r];
  }
  T(3,3) = R[0].one();
  return T;
}

PolyMatrix identity3(const Poly & zero)
{
  PolyMatrix I(zero,3,3);
  for( int i = 0; i < 3; ++i )
    I(i,i) = zero.one();
  return I;
}

PolyMatrix topLeftRotation(PolyMatrix & T)
{
  PolyMatrix R(T[0].zero(),3,3);
  for( int r = 0; r < 3; ++r )
    for( int c = 0; c < 3; ++c )
      R(r,c) = T(r,c);
  return R;
}

PolyMatrix rightTranslation(PolyMatrix & T)
{
  PolyMatrix t(T[0].zero(),3,1);
  for( int r = 0; r < 3; ++r )
    t[r] = T(r,3);
  return t;
}

PolyMatrix translationDesign(PolyMatrix & tB, size_t nu)
{
  // d_i(q)=t_Ai-R(q)t_Bi = t_Ai + D_i*r(q), where r(q) is row-major vec(R).
  PolyMatrix D(tB[0].zero(),3,9);
  for( int r = 0; r < 3; ++r )
    for( int c = 0; c < 3; ++c )
      D(r,ridx(r,c)) = tB[c].negation();
  return D;
}

PolyMatrix adjugate3(PolyMatrix & M)
{
  PolyMatrix A(M[0].zero(),3,3);

  Poly a = M(0,0);
  Poly b = M(0,1);
  Poly c = M(0,2);
  Poly d = M(1,0);
  Poly e = M(1,1);
  Poly f = M(1,2);
  Poly g = M(2,0);
  Poly h = M(2,1);
  Poly i = M(2,2);

  A(0,0) = e*i - f*h;
  A(0,1) = c*h - b*i;
  A(0,2) = b*f - c*e;

  A(1,0) = f*g - d*i;
  A(1,1) = a*i - c*g;
  A(1,2) = c*d - a*f;

  A(2,0) = d*h - e*g;
  A(2,1) = b*g - a*h;
  A(2,2) = a*e - b*d;

  return A;
}

Poly derivative(const Poly & f, int variableIndex)
{
  const size_t nu = f.leadingTerm().monomial().dimensions();
  Poly result = Poly::zeroSZ(nu);

  Poly::terms_t::iterator it = f.begin();
  while( it != f.end() )
  {
    Term t = it->clone();
    const unsigned int exponent = t.monomial().exponents()[variableIndex];

    if( exponent > 0 )
    {
      std::vector<unsigned int> newExponents = t.monomial().exponents();
      newExponents[variableIndex]--;

      Coefficient coeff_zp = t.coefficient().clone();
      t.setDominant(1);
      Coefficient coeff_sym = t.coefficient().clone();
      t.setDominant(0);

      result += Poly::constSZ(exponent,nu) *
          Poly(Term(coeff_sym,coeff_zp,Monomial(newExponents)));
    }

    ++it;
  }

  return result;
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

bool appendNonzero(vector<Poly> & coeffs, const Poly & value)
{
  if( value.isZero() )
    return false;
  coeffs.push_back(value.clone());
  return true;
}

bool buildReducedCoefficientSample(
    vector<PolyMatrix> & A,
    vector<PolyMatrix> & B,
    size_t nu,
    vector<Poly> & coeffs)
{
  coeffs.clear();

  const int numberPairs = static_cast<int>(A.size());
  const Poly zero = Poly::zeroZ(nu);
  const PolyMatrix I = identity3(zero);

  PolyMatrix H(zero,3,3);
  PolyMatrix gAlpha(zero,3,1);
  PolyMatrix G(zero,3,9);

  for( int i = 0; i < numberPairs; ++i )
  {
    PolyMatrix RA = topLeftRotation(A[i]);
    PolyMatrix tA = rightTranslation(A[i]);
    PolyMatrix tB = rightTranslation(B[i]);
    PolyMatrix C = RA - I;
    PolyMatrix D = translationDesign(tB,nu);

    H += C.transpose() * C;
    gAlpha += C.transpose() * tA;
    G += C.transpose() * D;
  }

  Poly detH = H.determinant();
  if( detH.isZero() )
    return false;

  PolyMatrix adjH = adjugate3(H);
  PolyMatrix yAlpha = (adjH * gAlpha).negation();
  PolyMatrix Y = (adjH * G).negation();

  for( int i = 0; i < numberPairs; ++i )
  {
    PolyMatrix RA = topLeftRotation(A[i]);
    PolyMatrix RB = topLeftRotation(B[i]);
    PolyMatrix tA = rightTranslation(A[i]);
    PolyMatrix tB = rightTranslation(B[i]);
    PolyMatrix C = RA - I;
    PolyMatrix D = translationDesign(tB,nu);

    // det(H) * vec(R_A R(q) - R(q) R_B), using only structurally non-zero
    // coefficients for each linear residual in vec(R(q)).
    for( int a = 0; a < 3; ++a )
    {
      for( int b = 0; b < 3; ++b )
      {
        for( int r = 0; r < 3; ++r )
        {
          for( int c = 0; c < 3; ++c )
          {
            if( c != b && r != a )
              continue;

            Poly value = zero;
            if( c == b )
              value += RA(a,r);
            if( r == a )
              value -= RB(c,b);

            if( !appendNonzero(coeffs,detH * value) )
              return false;
          }
        }
      }
    }

    // C_i*(-adj(H)g(q)) + det(H)*(t_Ai-R(q)t_Bi).
    PolyMatrix transAlpha = C * yAlpha + tA * detH;
    PolyMatrix transBeta = C * Y + D * detH;

    for( int k = 0; k < 3; ++k )
    {
      if( !appendNonzero(coeffs,transAlpha[k]) )
        return false;

      for( int m = 0; m < 9; ++m )
      {
        if( !appendNonzero(coeffs,transBeta(k,m)) )
          return false;
      }
    }
  }

  return true;
}

} // namespace

int main(int argc, char** argv)
{
  (void)argc;
  (void)argv;

  initGenerator();

  // The generated Groebner solver is quaternion-only. Translation is removed
  // analytically and reconstructed in the public C++ wrapper after q is found.
  const size_t nu = 4;
  const int numberPairs = POLYJAM_HAND_EYE_PAIRS;

  std::cout << "Generating the determinant-cleared quaternion-only hand-eye AX=XB "
            << "Groebner solver with " << numberPairs << " motion pairs." << std::endl;

  vector<PolyMatrix> AgtList;
  vector<PolyMatrix> BgtList;
  vector<Poly> coeffSample;

  bool generatedValidSample = false;
  for( int attempt = 0; attempt < 512 && !generatedValidSample; ++attempt )
  {
    AgtList.clear();
    BgtList.clear();
    AgtList.reserve(numberPairs);
    BgtList.reserve(numberPairs);

    // Ground-truth X = [R(q_X), t_X] in Zp. The synthetic A_i/B_i coefficients
    // are exactly consistent with this X, which keeps the ideal non-empty in
    // the finite-field analysis instance.
    PolyMatrix qXgt = randomUnitQuaternionZ(nu);
    PolyMatrix RXgt = quaternionRotationZ(qXgt[0],qXgt[1],qXgt[2],qXgt[3]);
    PolyMatrix tXgt = randomTranslationZ(nu);

    for( int i = 0; i < numberPairs; ++i )
    {
      PolyMatrix qBgt = randomUnitQuaternionZ(nu);
      PolyMatrix RBgt = quaternionRotationZ(qBgt[0],qBgt[1],qBgt[2],qBgt[3]);
      PolyMatrix tBgt = randomTranslationZ(nu);

      // AX = XB implies:
      //   R_A = R_X R_B R_X^T
      //   t_A = R_X t_B + t_X - R_A t_X
      PolyMatrix RAgt = RXgt * RBgt * RXgt.transpose();
      PolyMatrix tAgt = RXgt * tBgt + tXgt - RAgt * tXgt;

      AgtList.push_back(makeTransform(RAgt,tAgt,nu));
      BgtList.push_back(makeTransform(RBgt,tBgt,nu));
    }

    generatedValidSample = buildReducedCoefficientSample(
        AgtList,BgtList,nu,coeffSample);
  }

  if( !generatedValidSample )
  {
    std::cerr << "Unable to generate a finite-field hand-eye instance with "
              << "non-zero determinant-cleared reduced-objective coefficients."
              << std::endl;
    return 1;
  }

  std::cout << "Using " << coeffSample.size()
            << " scalar coefficients for the reduced objective."
            << std::endl;

  // Unknown rotation q = (qw,qx,qy,qz). The public wrapper returns the final
  // [qw,qx,qy,qz,tx,ty,tz]^T solution after sequentially reconstructing t.
  Poly qw = Poly::uSZ(1,nu);
  Poly qx = Poly::uSZ(2,nu);
  Poly qy = Poly::uSZ(3,nu);
  Poly qz = Poly::uSZ(4,nu);
  PolyMatrix q(Poly::zeroSZ(nu),4,1);
  q[0] = qw;
  q[1] = qx;
  q[2] = qy;
  q[3] = qz;

  PolyMatrix R = quaternionRotation(qw,qx,qy,qz);
  PolyMatrix r(Poly::zeroSZ(nu),9,1);
  for( int row = 0; row < 3; ++row )
    for( int col = 0; col < 3; ++col )
      r[ridx(row,col)] = R(row,col);

  // The original scalar Euclidean objective is
  //
  //   E(q,t) = sum_i || R_Ai R(q) - R(q) R_Bi ||_F^2
  //          + sum_i || (R_Ai-I)t + t_Ai - R(q)t_Bi ||_2^2.
  //
  // The translation Jacobian gives Ht+g(q)=0.  The public wrapper computes the
  // same determinant-cleared residual coefficients represented here by coeffs.
  // This avoids symbolic division and keeps the generated template compact:
  //
  //   Ebar(q) = sum_k residual_k(q)^2,
  //
  // where the residuals are linear in vec(R(q)) after substituting
  // det(H)t(q)=-adj(H)g(q).
  Poly reducedObjective = Poly::zeroSZ(nu);
  size_t coeffIndex = 0;

  for( int pair = 0; pair < numberPairs; ++pair )
  {
    for( int a = 0; a < 3; ++a )
    {
      for( int b = 0; b < 3; ++b )
      {
        Poly residual = Poly::zeroSZ(nu);
        for( int row = 0; row < 3; ++row )
        {
          for( int col = 0; col < 3; ++col )
          {
            if( col != b && row != a )
              continue;
            residual += liftCoefficient(coeffSample[coeffIndex],coeffIndex,nu) *
                        r[ridx(row,col)];
            ++coeffIndex;
          }
        }
        reducedObjective += residual * residual;
      }
    }

    for( int k = 0; k < 3; ++k )
    {
      Poly residual = liftCoefficient(coeffSample[coeffIndex],coeffIndex,nu);
      ++coeffIndex;

      for( int m = 0; m < 9; ++m )
      {
        residual += liftCoefficient(coeffSample[coeffIndex],coeffIndex,nu) * r[m];
        ++coeffIndex;
      }

      reducedObjective += residual * residual;
    }
  }

  if( coeffIndex != coeffSample.size() )
  {
    std::cerr << "Internal coefficient-count mismatch." << std::endl;
    return 2;
  }

  Poly unitResidual = qw*qw + qx*qx + qy*qy + qz*qz - Poly::oneSZ(nu);

  Poly grad[4] = {
    derivative(reducedObjective,0),
    derivative(reducedObjective,1),
    derivative(reducedObjective,2),
    derivative(reducedObjective,3)
  };

  std::list<Poly*> eqs;
  eqs.push_back(new Poly(unitResidual));

  // Criticality on S^3: grad(Ebar) must be parallel to q. Eliminating the
  // Lagrange multiplier produces the six 2x2 minors q_i*grad_j-q_j*grad_i.
  for( int i = 0; i < 4; ++i )
  {
    for( int j = i + 1; j < 4; ++j )
      eqs.push_back(new Poly(q[i] * grad[j] - q[j] * grad[i]));
  }

  std::cout << "Constructed a scalar objective, eliminated translation with "
            << "the translation Jacobian, and formed " << eqs.size()
            << " quaternion-only critical equations in " << nu
            << " unknowns." << std::endl;

  const string parameters("const std::vector<double> & coeffs");

  execGenerator(eqs,string("handeye_axxb_q"),parameters,false);

  for( std::list<Poly*>::iterator it = eqs.begin(); it != eqs.end(); ++it )
    delete *it;

  return 0;
}
