/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** ****************************************************************** */

// Written: Javad Esmaeelpour
// Created: 08/2026

#include <ExtrapolationPredictor.h>
#include <AnalysisModel.h>
#include <IncrementalIntegrator.h>
#include <Vector.h>
#include <math.h>

ExtrapolationPredictor::ExtrapolationPredictor(double g, int minS)
 :gain(g), minSteps(minS), numObserved(0),
  dUprev(0), dUprev2(0), lastConfidence(0.0)
{

}

ExtrapolationPredictor::~ExtrapolationPredictor()
{
    if (dUprev  != 0) delete dUprev;
    if (dUprev2 != 0) delete dUprev2;
}

void
ExtrapolationPredictor::reset(void)
{
    if (dUprev  != 0) { delete dUprev;  dUprev  = 0; }
    if (dUprev2 != 0) { delete dUprev2; dUprev2 = 0; }
    numObserved    = 0;
    lastConfidence = 0.0;
}

bool
ExtrapolationPredictor::predict(const AnalysisModel &,
                                const IncrementalIntegrator &,
                                Vector &dU_guess)
{
    lastConfidence = 0.0;

    if (numObserved < minSteps || dUprev == 0 || dUprev2 == 0)
        return false;
    if (dUprev->Size() != dU_guess.Size() || dUprev2->Size() != dU_guess.Size())
        return false;

    double n1 = dUprev->Norm();
    double n2 = dUprev2->Norm();

    // Both increments essentially zero: nothing is happening and the
    // integrator's predicted point is already right. Declining here
    // matters -- it avoids paying the accept/reject overhead on steps
    // with nothing to gain.
    if (n1 < 1.0e-14 || n2 < 1.0e-14)
        return false;

    // dU_guess = gain * (dU_{n-1} - dU_{n-2})
    dU_guess = *dUprev;
    dU_guess.addVector(1.0, *dUprev2, -1.0);
    dU_guess *= gain;

    // Two independent smoothness signals, combined multiplicatively so
    // that either one going bad kills the guess.
    //
    //   alignment  cosine between consecutive increments. Near 1 the
    //              response tracks a consistent direction; a yield
    //              event, unloading reversal or gap closure swings it.
    //   steadiness ratio of consecutive magnitudes, folded to (0,1].
    //              Near 1 means a stable step-to-step rate.
    double dot = (*dUprev) ^ (*dUprev2);
    double alignment = dot / (n1 * n2);
    if (alignment < 0.0)
        alignment = 0.0;                 // direction reversed -> no trust

    double steadiness = (n1 < n2) ? (n1 / n2) : (n2 / n1);

    lastConfidence = alignment * steadiness;

    return true;
}

void
ExtrapolationPredictor::observe(const AnalysisModel &,
                                const IncrementalIntegrator &,
                                const Vector &dU_converged)
{
    // NOTE: assumes a constant step size. Under variable dt, or an
    // adaptive integrator that subdivides, consecutive increments are
    // not directly comparable and this extrapolation is biased. The fix
    // is to scale each stored increment by its own dt before
    // differencing. Left undone deliberately: benchmark at constant dt
    // first, then decide whether variable dt is worth supporting.

    if (dUprev2 != 0) delete dUprev2;
    dUprev2 = dUprev;
    dUprev  = new Vector(dU_converged);

    numObserved++;
}
