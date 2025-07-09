//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//----------------------------------------------------------------------------//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Please cite the following resource in any derivative works:
//
// [1] Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//

//
// Description: SouzaFrameTransf implements a Corotational
// transformation for a spatial frame element following the formulation
// by Crisfield (1990).
//
// Written: Claudio Perez
// Created: 05/2024
//
// Adapted from: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
//
// [2] Crisfield, M.A. (1990) "A consistent co-rotational formulation for
//     non-linear, three-dimensional, beam-elements", Computer Methods in Applied
//     Mechanics and Engineering, 81(2), pp. 131–150. Available at:
//     https://doi.org/10.1016/0045-7825(90)90106-V.
//
#include <math.h>
#include <Node.h>
#include <Logging.h>
#include "SouzaFrameTransf.h"

#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <GroupSO3.h>
#include "Isometry/CrisfieldIsometry.h"

namespace OpenSees {

template <int nn, int ndf>
SouzaFrameTransf<nn,ndf>::SouzaFrameTransf(int tag, const Vector3D &vz,
                                            const std::array<Vector3D, nn> *offset,
                                            int offset_flags)
  : FrameTransform<nn,ndf>(tag),
    vz(vz),
    offsets{nullptr},
    L(0), Ln(0),
    nodeIInitialDisp(0), nodeJInitialDisp(0),
    initialDispChecked(false),
    crs{vz}
{
  alphaI.zero();
  alphaJ.zero();

  // Rigid joint offsets
  if (offset != nullptr) {
    offsets = new std::array<Vector3D, nn>{};
    *offsets = *offset;
  }
}


template <int nn, int ndf>
SouzaFrameTransf<nn,ndf>::~SouzaFrameTransf()
{
  if (offsets != nullptr)
    delete offsets;
}


template <int nn, int ndf>
double
SouzaFrameTransf<nn,ndf>::getInitialLength()
{
  return L;
}


template <int nn, int ndf>
double
SouzaFrameTransf<nn,ndf>::getDeformedLength()
{
  return Ln;
}


template <int nn, int ndf>
FrameTransform<nn,ndf> *
SouzaFrameTransf<nn,ndf>::getCopy() const
{

  SouzaFrameTransf *theCopy =
    new SouzaFrameTransf<nn,ndf>(this->getTag(), vz, offsets);

  theCopy->nodes[0]  = nodes[0];
  theCopy->nodes[1]  = nodes[1];
  theCopy->xAxis     = xAxis;
  theCopy->L         = L;
  theCopy->Ln        = Ln;
  theCopy->R0        = R0;
  theCopy->Q_pres[0] = Q_pres[0];
  theCopy->Q_pres[1] = Q_pres[1];
  theCopy->Q_past[0] = Q_past[0];
  theCopy->Q_past[1] = Q_past[1];
  theCopy->ul = ul;
  theCopy->ulcommit = ulcommit;
  return theCopy;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::revertToStart()
{
  ul.zero();
  Q_pres[0] = Versor::from_matrix(R0);
  for (int i=1; i<nn; i++)
    Q_pres[i] = Q_pres[0];

  alphaI.zero();
  alphaJ.zero();

  this->update();
  return 0;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::commit()
{
  ulcommit  = ul;
  Q_past[0] = Q_pres[0];
  Q_past[1] = Q_pres[1];
  return 0;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::revertToLastCommit()
{
  // determine global displacement increments from last iteration
  const Vector &dispI = nodes[0]->getTrialDisp();
  const Vector &dispJ = nodes[1]->getTrialDisp();

  for (int k = 0; k < 3; k++) {
    alphaI(k) =  dispI(k+3);
    alphaJ(k) =  dispJ(k+3);
  }

  ul        = ulcommit;
  Q_pres[0] = Q_past[0];
  Q_pres[1] = Q_past[1];

  this->update();

  return 0;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::initialize(std::array<Node*, nn>& new_nodes)
{
  for (int i=0; i<nn; i++) {
    nodes[i] = new_nodes[i];
    if (nodes[i] == nullptr) {
      opserr << "invalid pointers to the element nodes\n";
      return -1;
    }
  }

  crs.initialize(new_nodes);

  dX  = nodes[nn-1]->getCrds();
  dX -= nodes[   0]->getCrds();

  // Add initial displacements at nodes
  if (initialDispChecked == false) {
    const Vector &nodeIDisp = nodes[0]->getDisp();
    const Vector &nodeJDisp = nodes[1]->getDisp();
    for (int i = 0; i<6; i++)
      if (nodeIDisp[i] != 0.0) {
        nodeIInitialDisp = new double [6];
        for (int j = 0; j<6; j++)
          nodeIInitialDisp[j] = nodeIDisp[j];
        i = 6;
      }

    for (int j = 0; j<6; j++)
      if (nodeJDisp[j] != 0.0) {
        nodeJInitialDisp = new double [6];
        for (int i = 0; i<6; i++)
          nodeJInitialDisp[i] = nodeJDisp[i];
        j = 6;
      }
    initialDispChecked = true;
  }

  //
  // Length and Orientation
  //
  Vector3D dx;

  dx = nodes[nn-1]->getCrds() - nodes[0]->getCrds();

  L = dx.norm();

  if (L == 0.0) {
    opserr << "\nSouzaFrameTransf::computeElemtLengthAndOrien: 0 length\n";
    return -2;
  }

  //
  // Set rotation matrix
  //
  int error = FrameTransform<nn,ndf>::Orient(dx, vz, R0);
  if (error)
    return error;

  // Compute initial pseudo-vectors for nodal triads
  Q_pres[0] = Q_pres[1] = Versor::from_matrix(R0);

  ul.zero();
  ulpr.zero();

  for (int i=0; i<nn; i++)
    vr[i].zero();

  this->commit();

  return 0;
}

template <int nn, int ndf>
VectorND<nn*ndf>
SouzaFrameTransf<nn,ndf>::getStateVariation()
{
  return ul - ulpr;
}


template <int nn, int ndf>
Vector3D
SouzaFrameTransf<nn,ndf>::getNodePosition(int tag)
{
  Vector3D u;
  for (int i=0; i<3; i++)
    u[i] = ul[tag*ndf+i];
  return u;
}


template <int nn, int ndf>
Vector3D
SouzaFrameTransf<nn,ndf>::getNodeRotationLogarithm(int tag)
{
  return vr[tag];
}


//
// Set RI,RJ,Rbar, Ln, e and ul
//
template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::update()
{
  // determine global displacement increments from last iteration

  const Vector& dispI = nodes[   0]->getTrialDisp();
  const Vector& dispJ = nodes[nn-1]->getTrialDisp();

  //
  // Update state
  //

  // 1.1 Relative translation
  Vector3D dx = dX;// dx = dX + dJI;
  {
    // Relative translation
    for (int k = 0; k < 3; k++)
      dx[k] += dispJ(k) - dispI(k);


    // Calculate the deformed length
    Ln = dx.norm();

    if (Ln == 0.0) {
      opserr << "\nSouzaFrameTransf: deformed length is 0.0\n";
      return -2;
    }
  }

  // 1.2 Rotational displacement increments
  {
    Vector3D dAlphaI, dAlphaJ;

    for (int k = 0; k < 3; k++) {
      dAlphaI[k] =  dispI(k+3) - alphaI[k];
      alphaI[k]  =  dispI(k+3);
      dAlphaJ[k] =  dispJ(k+3) - alphaJ[k];
      alphaJ[k]  =  dispJ(k+3);
    }

    // Update the nodal rotations; Note the Hamilton product is assumed!
    if (dAlphaI.norm() != 0)
      Q_pres[0] = Versor::from_vector(dAlphaI)*Q_pres[0];
    if (dAlphaJ.norm() != 0)
      Q_pres[1] = Versor::from_vector(dAlphaJ)*Q_pres[1];
  }

  //
  // 2) Form transformation
  //

  crs.update(MatrixFromVersor(Q_pres[0]), 
             MatrixFromVersor(Q_pres[1]), 
             dx, nodes);

  //
  // 3) Local deformations
  //

  // Save previous state
  ulpr = ul;

  // Rotations
  {
    Matrix3D e = crs.getRotation();
    vr[0] = LogC90(e^MatrixFromVersor(Q_pres[0]));
    for (int i=0; i<3; i++)
      ul[imx+i] = vr[0][i];

    vr[1] = LogC90(e^MatrixFromVersor(Q_pres[1]));
    for (int i=0; i<3; i++)
      ul[jmx+i] = vr[1][i];
  }

  // Axial
  ul(inx) = 0;
  ul(jnx) = Ln - L;
  // ul(jnx) = (dX + (dx-dX)*0.5).dot(dx-dX)*2.0/(Ln+L);

  // Form the transformation tangent
  T = crs.compute_tangent(ul);
  return 0;
}


template <int nn, int ndf>
// inline VectorND<nn*ndf>
int
SouzaFrameTransf<nn,ndf>::push(VectorND<nn*ndf>&pl, Operation)
{
  // return T^pl;
  VectorND<nn*ndf> pg{};
  for (int a = 0; a<nn; a++) {
    VectorND<6> pa {pl(a*ndf+0), pl(a*ndf+1), pl(a*ndf+2), 
                    pl(a*ndf+3), pl(a*ndf+4), pl(a*ndf+5)};

    for (int b = 0; b<2; b++) {
      VectorND<6> pab{};
      for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
          pab[j] += T(a*6 + i, b*6+j) * pa[i];
      // VectorND<6> pab = pushResponse(pa, a, b);
      pg.assemble(b*6, pab, 1.0);
    }
  }
  pl = pg; // TODO: optimize
  return 0;
}


// do 
//    K = ag'*Km*ag + Kp
//
template <int nn, int ndf>
// MatrixND<nn*ndf,nn*ndf>
int
SouzaFrameTransf<nn,ndf>::push(MatrixND<nn*ndf,nn*ndf>& kl, 
                               const VectorND<nn*ndf>& pl,
                               Operation op)
{    
  static MatrixND<12,12> K;
  K.addMatrixTripleProduct(0.0, T, kl, 1.0);

  // Add geometric part kg
  this->addTangent(K, pl, ul);

  kl = K; // TODO: optimize
  return 0;
}


//
// Add geometric part of the transformation tangent
//
//  kg += T'*kl*T + ks1 + T * diag(m.*tan(thetal))*T' + ...
//         m(4)*(ks2r2t3_u3 + ks2r3u2_t2) + ...
//         m(2)*ks2r2t1 + m(3)*ks2r3t1 + ...
//         m(5)*ks2r2u1 + m(6)*ks2r3u1 + ...
//         ks3 + ks3' + ks4 + ks5;
template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::addTangent(MatrixND<12,12>& kg,
                                     const VectorND<12>& pl,
                                     const VectorND<12>& ul)
{

  crs.addTangent(kg, pl, ul);

  //
  //  T' * diag (M .* tan(thetal))*T
  //
  for (int node=0; node<2; node++) {
    for (int k = 0; k < 3; k++) {
      const double factor =  pl[(node ? jmx : imx) + k] // pl[6*node+3+k]
                          * std::tan(ul[(node ? jmx : imx) + k]);


      for (int i = 0; i < 12; i++) {
        const double Tki = T((node ? jmx : imx) + k,i);
        for (int j = 0; j < 12; j++) {
          kg(i,j) += Tki * factor * T((node ? jmx : imx) + k, j);
        }
      }
    }
  }

  return 0;
}


template <int nn, int ndf>
int
SouzaFrameTransf<nn,ndf>::getLocalAxes(Vector3D &e1, Vector3D &e2, Vector3D &e3) const
{
  for (int i = 0; i < 3; i++) {
    e1[i] = R0(i,0);
    e2[i] = R0(i,1);
    e3[i] = R0(i,2);
  }
  return 0;
}


template <int nn, int ndf>
double
SouzaFrameTransf<nn,ndf>::getLengthGrad()
{
  const int di = nodes[0]->getCrdsSensitivity();
  const int dj = nodes[1]->getCrdsSensitivity();

  Vector3D dxi{0.0};
  Vector3D dxj{0.0};

  if (di != 0)
    dxi(di-1) = 1.0;
  if (dj != 0)
    dxj(dj-1) = 1.0;

  return 1.0/L * dX.dot(dxj - dxi);
}

template <int nn, int ndf>
const Vector &
SouzaFrameTransf<nn,ndf>::getBasicDisplTotalGrad(int gradNumber)
{
  opserr << "WARNING CrdTransf::getBasicDisplTotalGrad() - this method "
         << " should not be called.\n";

  static Vector dummy(1);
  return dummy;
}

template <int nn, int ndf>
const Vector &
SouzaFrameTransf<nn,ndf>::getBasicDisplFixedGrad()
{
  opserr << "ERROR CrdTransf::getBasicDisplFixedGrad() - has not been"
          << " implemented yet for the chosen transformation\n.";

  static Vector dummy(1);
  return dummy;
}


template <int nn, int ndf>
const Vector &
SouzaFrameTransf<nn,ndf>::getGlobalResistingForceShapeSensitivity(const Vector &pb,
						   const Vector &p0,
						   int gradNumber)
{
  opserr << "ERROR CrdTransf::getGlobalResistingForceSensitivity() - has not been"
      << " implemented yet for the chosen transformation." << endln;
  
  static Vector dummy(1);
  return dummy;
}

template <int nn, int ndf>
void
SouzaFrameTransf<nn,ndf>::Print(OPS_Stream &s, int flag)
{

  if (flag == OPS_PRINT_CURRENTSTATE) {
      s << "\nFrameTransform: " << this->getTag() << " Type: SouzaFrameTransf";
      s << "\tvxz: " << Vector(vz);
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"SouzaFrameTransf\"" << ", ";
    s << "\"vecxz\": [" << vz(0) << ", " << vz(1) << ", " << vz(2) << "]";

    if (offsets != nullptr) {
      s << ", \"offsets\": [";
      for (int i=0; i<nn; i++) {
        s << "["
          << (*offsets)[i][0] << ", " 
          << (*offsets)[i][1] << ", "
          << (*offsets)[i][2] << "]";
        if (i < nn-1)
          s << ", ";
      }
      s << "]";
    }
    s << "}";
    return;
  }
}

} // namespace OpenSees