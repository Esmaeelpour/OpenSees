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
//
// Description: the deterministic baseline NewtonPredictor for
// WarmStartNewton. No external dependency, no training step, identical
// output on every run.
//
// Method: from the last two converged increments dU_{n-1}, dU_{n-2},
// predict the correction to the integrator's own predicted point as
//
//     dU_guess = gain * (dU_{n-1} - dU_{n-2})
//
// i.e. continue the observed trend in the increment. Doing nothing
// (order zero) is what OpenSees does today.
//
// Confidence comes from how steady the recent history is. Nearly
// parallel, similarly sized consecutive increments mean the response is
// on a smooth branch and extrapolation is reliable; a sharp change in
// direction or magnitude means a state event just occurred, so
// confidence drops and the trust gate rejects the guess.

#ifndef ExtrapolationPredictor_h
#define ExtrapolationPredictor_h

#include <WarmStartNewton.h>
#include <Vector.h>

class ExtrapolationPredictor : public NewtonPredictor
{
public:
    // gain     multiplier on the extrapolated trend; 1.0 is the natural
    //          linear extrapolation, below 1.0 is damped. Tune only
    //          from benchmark evidence.
    // minSteps converged steps required before predicting at all.
    ExtrapolationPredictor(double gain = 1.0, int minSteps = 3);
    ~ExtrapolationPredictor();

    bool predict(const AnalysisModel &theModel,
                 const IncrementalIntegrator &theIntegrator,
                 Vector &dU_guess);

    void observe(const AnalysisModel &theModel,
                 const IncrementalIntegrator &theIntegrator,
                 const Vector &dU_converged);

    double confidence(void) const { return lastConfidence; }

    void reset(void);

    const char *getType(void) const { return "Extrapolation"; }

private:
    double gain;
    int minSteps;
    int numObserved;

    Vector *dUprev;    // dU_{n-1}
    Vector *dUprev2;   // dU_{n-2}

    double lastConfidence;
};

#endif
