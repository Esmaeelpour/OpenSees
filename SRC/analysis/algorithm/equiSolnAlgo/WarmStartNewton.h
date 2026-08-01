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
// Description: WarmStartNewton is a Newton-Raphson algorithm that, before
// entering the iteration loop, asks a pluggable NewtonPredictor for a
// better starting point than the integrator's own predicted step.
//
// It does NOT replace Newton's local convergence guarantee. The
// formTangent/solve/update/formUnbalance loop is that of NewtonRaphson
// and remains the correctness backbone; convergence is still decided
// entirely by the user's ConvergenceTest, so the converged answer is
// unchanged and only the path to it is (hopefully) cheaper.
//
// Acceptance is never-worse: a guess is applied only if it strictly
// reduces the residual norm relative to the integrator's predicted
// point, and is otherwise undone exactly. Cost of that guard is one
// extra formUnbalance() and up to two extra update() calls per step, so
// a predictor with a poor hit rate can be a net loss -- check the
// diagnostics counters before trusting one on a production run.

#ifndef WarmStartNewton_h
#define WarmStartNewton_h

#include <EquiSolnAlgo.h>
#include <Vector.h>

class AnalysisModel;
class IncrementalIntegrator;
class LinearSOE;

// ---------------------------------------------------------------------
// NewtonPredictor -- pluggable warm-start interface.
//
// Deliberately framework-agnostic: the reference implementation
// (ExtrapolationPredictor) is plain polynomial extrapolation with no
// external dependency, but the same interface admits a POD/PGD reduced
// basis or a learned model without OpenSees core taking on any of it.
// ---------------------------------------------------------------------
class NewtonPredictor
{
public:
    virtual ~NewtonPredictor() {}

    // Predicted displacement increment for the step about to be solved,
    // measured RELATIVE to where the integrator's newStep() predictor
    // has already placed the model. Return false to decline.
    //
    // Model and integrator are const so the "identical when declined"
    // guarantee is enforced by the compiler, not by comment.
    virtual bool predict(const AnalysisModel &theModel,
                         const IncrementalIntegrator &theIntegrator,
                         Vector &dU_guess) = 0;

    // Called once per converged step so an online predictor can update.
    virtual void observe(const AnalysisModel &theModel,
                         const IncrementalIntegrator &theIntegrator,
                         const Vector &dU_converged) {}

    // Confidence in [0,1] for the last predict(). Default is 0.0, NOT
    // 0.5: a predictor that does not estimate its own confidence must
    // not be trusted by default. Silence means "do not use me".
    virtual double confidence(void) const { return 0.0; }

    // Discard accumulated history; called on domainChanged().
    virtual void reset(void) {}

    virtual const char *getType(void) const = 0;
};

class WarmStartNewton: public EquiSolnAlgo
{
public:
    WarmStartNewton(int theTangentToUse = CURRENT_TANGENT,
                    double iFactor = 0.0,
                    double cFactor = 1.0,
                    double trustRadius = 0.25,
                    double minConfidence = 0.0,
                    double lookaheadFactor = 1.0,
                    NewtonPredictor *thePredictor = 0);
    WarmStartNewton();
    ~WarmStartNewton();

    // Takes ownership: the algorithm deletes the predictor.
    void setPredictor(NewtonPredictor *thePredictor);

    int solveCurrentStep(void);
    int setConvergenceTest(ConvergenceTest *theNewTest);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    int getNumIterations(void);

    // Diagnostics. The whole case for this algorithm rests on these:
    // report ITERATIONS SAVED, not merely how often the predictor fired.
    int getNumPredictedSteps(void) const { return numPredictedSteps; }
    int getNumFallbackSteps(void) const { return numFallbackSteps; }
    int getNumRejectedGuesses(void) const { return numRejectedGuesses; }
    double getAvgIterationsPredicted(void) const;
    double getAvgIterationsFallback(void) const;

protected:
    int domainChanged(void);

private:
    // Reference scale for the magnitude gate: the norm of the last
    // converged step increment. Neither AnalysisModel nor
    // IncrementalIntegrator exposes the integrator's own predicted
    // increment, so there is nothing else generic to scale against --
    // and this is the better choice anyway, being defined identically
    // under static and dynamic integrators.
    double referenceScale(void) const;

    int tangent;
    double iFactor, cFactor;
    double trustRadius;
    double minConfidence;

    // A single trial-point residual check cannot tell whether the guess
    // actually reduced the number of iterations still needed -- a guess
    // can look locally better yet leave the state somewhere Newton has
    // to spend MORE total iterations correcting (overshoot, wrong-
    // direction curvature). Guard: after the guess passes the residual
    // check, form the tangent and solve for the correction AT the
    // guessed point -- this is the mandatory first Newton iteration
    // regardless, not wasted work. Only truly commit to the guess if
    // that correction is no larger than the guess itself
    // (||dX|| <= lookaheadFactor * ||dUguess||): if fixing up after the
    // guess needs a correction as big as the guess, the guess bought
    // nothing and is more likely to have overshot or moved off in a bad
    // direction than to be genuine progress.
    double lookaheadFactor;

    int numIterations;
    int numPredictedSteps;
    int numFallbackSteps;
    int numRejectedGuesses;
    long totalIterationsPredicted;
    long totalIterationsFallback;

    // The displacement change this algorithm applies over a step, i.e.
    // the correction on top of whatever the integrator's newStep()
    // predictor already did. Accumulated from the solved increments
    // (and any accepted guess) because no generic accessor for it
    // exists. This -- not the total nodal displacement change -- is
    // what the predictor learns and predicts, which is exactly right:
    // warm-starting is about anticipating the correction newStep()
    // leaves behind.
    Vector *stepIncr;
    Vector *lastConvergedIncr;
    Vector *guessWork;

    NewtonPredictor *thePredictor;   // owned; null -> pure NewtonRaphson
};

#endif
