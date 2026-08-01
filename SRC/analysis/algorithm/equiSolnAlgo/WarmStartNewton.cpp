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
// Description: implementation of WarmStartNewton. The iteration loop is
// deliberately kept structurally identical to NewtonRaphson.cpp so that
// a diff of the two shows the warm-start block as the only addition.

#include <WarmStartNewton.h>
#include <ExtrapolationPredictor.h>
#include <AnalysisModel.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <ConvergenceTest.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>

#include <elementAPI.h>
#include <string.h>
#include <math.h>

void *OPS_WarmStartNewton()
{
    int formTangent = CURRENT_TANGENT;
    double iFactor = 0.0, cFactor = 1.0;
    double trustRadius = 0.25;
    double minConfidence = 0.0;
    double lookaheadFactor = 1.0;

    // Predictor parameters. Extrapolation is the default and the only
    // predictor built in -- it is deterministic and dependency-free, and
    // is the baseline any other predictor has to beat.
    double gain = 1.0;
    int minSteps = 3;
    bool wantPredictor = true;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *type = OPS_GetString();

        if ((strcmp(type, "-initial") == 0) || (strcmp(type, "-Initial") == 0) ||
            (strcmp(type, "-initialTangent") == 0)) {
            formTangent = INITIAL_TANGENT;

        } else if ((strcmp(type, "-currentTangent") == 0) ||
                   (strcmp(type, "-current") == 0)) {
            formTangent = CURRENT_TANGENT;

        } else if (strcmp(type, "-none") == 0) {
            // No predictor: reduces exactly to NewtonRaphson. Useful as
            // an in-build control when benchmarking, since it isolates
            // the predictor's effect from any other difference.
            wantPredictor = false;

        } else if (strcmp(type, "-predictor") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING WarmStartNewton -- -predictor needs a type\n";
                return 0;
            }
            const char *pType = OPS_GetString();
            if (strcmp(pType, "Extrapolation") == 0 ||
                strcmp(pType, "extrapolation") == 0) {
                wantPredictor = true;
            } else if (strcmp(pType, "none") == 0 || strcmp(pType, "None") == 0) {
                wantPredictor = false;
            } else {
                opserr << "WARNING WarmStartNewton -- unknown predictor '"
                       << pType << "'\n";
                return 0;
            }

        } else if (strcmp(type, "-gain") == 0) {
            int numData = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 ||
                OPS_GetDoubleInput(&numData, &gain) < 0) {
                opserr << "WARNING WarmStartNewton -- error reading gain\n";
                return 0;
            }

        } else if (strcmp(type, "-minSteps") == 0) {
            int numData = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 ||
                OPS_GetIntInput(&numData, &minSteps) < 0) {
                opserr << "WARNING WarmStartNewton -- error reading minSteps\n";
                return 0;
            }

        } else if (strcmp(type, "-trustRadius") == 0) {
            int numData = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 ||
                OPS_GetDoubleInput(&numData, &trustRadius) < 0) {
                opserr << "WARNING WarmStartNewton -- error reading trustRadius\n";
                return 0;
            }

        } else if (strcmp(type, "-minConfidence") == 0) {
            int numData = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 ||
                OPS_GetDoubleInput(&numData, &minConfidence) < 0) {
                opserr << "WARNING WarmStartNewton -- error reading minConfidence\n";
                return 0;
            }

        } else if (strcmp(type, "-lookaheadFactor") == 0) {
            int numData = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 ||
                OPS_GetDoubleInput(&numData, &lookaheadFactor) < 0) {
                opserr << "WARNING WarmStartNewton -- error reading lookaheadFactor\n";
                return 0;
            }
        }
    }

    NewtonPredictor *pred = 0;
    if (wantPredictor)
        pred = new ExtrapolationPredictor(gain, minSteps);

    return new WarmStartNewton(formTangent, iFactor, cFactor,
                               trustRadius, minConfidence, lookaheadFactor, pred);
}

WarmStartNewton::WarmStartNewton(int theTangentToUse, double iFact, double cFact,
                                 double trustRad, double minConf,
                                 double lookaheadFact,
                                 NewtonPredictor *pred)
 :EquiSolnAlgo(EquiALGORITHM_TAGS_WarmStartNewton),
  tangent(theTangentToUse), iFactor(iFact), cFactor(cFact),
  trustRadius(trustRad), minConfidence(minConf),
  lookaheadFactor(lookaheadFact),
  numIterations(0), numPredictedSteps(0), numFallbackSteps(0),
  numRejectedGuesses(0),
  totalIterationsPredicted(0), totalIterationsFallback(0),
  stepIncr(0), lastConvergedIncr(0), guessWork(0), thePredictor(pred)
{

}

WarmStartNewton::WarmStartNewton()
 :EquiSolnAlgo(EquiALGORITHM_TAGS_WarmStartNewton),
  tangent(CURRENT_TANGENT), iFactor(0.0), cFactor(1.0),
  trustRadius(0.25), minConfidence(0.0),
  lookaheadFactor(1.0),
  numIterations(0), numPredictedSteps(0), numFallbackSteps(0),
  numRejectedGuesses(0),
  totalIterationsPredicted(0), totalIterationsFallback(0),
  stepIncr(0), lastConvergedIncr(0), guessWork(0), thePredictor(0)
{

}

WarmStartNewton::~WarmStartNewton()
{
    if (thePredictor != 0)      delete thePredictor;
    if (stepIncr != 0)          delete stepIncr;
    if (lastConvergedIncr != 0) delete lastConvergedIncr;
    if (guessWork != 0)         delete guessWork;
}

void
WarmStartNewton::setPredictor(NewtonPredictor *pred)
{
    if (thePredictor != 0)
        delete thePredictor;
    thePredictor = pred;
}

int
WarmStartNewton::setConvergenceTest(ConvergenceTest *newTest)
{
    theTest = newTest;
    return 0;
}

int
WarmStartNewton::getNumIterations(void)
{
    return numIterations;
}

double
WarmStartNewton::getAvgIterationsPredicted(void) const
{
    if (numPredictedSteps == 0) return 0.0;
    return (double)totalIterationsPredicted / (double)numPredictedSteps;
}

double
WarmStartNewton::getAvgIterationsFallback(void) const
{
    if (numFallbackSteps == 0) return 0.0;
    return (double)totalIterationsFallback / (double)numFallbackSteps;
}

double
WarmStartNewton::referenceScale(void) const
{
    if (lastConvergedIncr != 0) {
        double s = lastConvergedIncr->Norm();
        if (s > 1.0e-14)
            return s;
    }
    // No usable scale yet: skip the magnitude gate and rely on the
    // residual accept/reject, which is the stronger check anyway.
    return 0.0;
}

int
WarmStartNewton::solveCurrentStep(void)
{
    AnalysisModel         *theAnaModel   = this->getAnalysisModelPtr();
    IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
    LinearSOE             *theSOE        = this->getLinearSOEptr();

    if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0)
        || (theTest == 0)) {
        opserr << "WARNING WarmStartNewton::solveCurrentStep() - setLinks() has";
        opserr << " not been called - or no ConvergenceTest has been set\n";
        return -5;
    }

    if (theIntegrator->formUnbalance() < 0) {
        opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
        opserr << "the Integrator failed in formUnbalance()\n";
        return -2;
    }

    // ------------------------------------------------------------------
    // Warm start. Skipped entirely when no predictor is registered, in
    // which case what follows is NewtonRaphson.
    // ------------------------------------------------------------------
    bool usedPrediction = false;

    int n = theSOE->getNumEqn();

    if (stepIncr == 0 || stepIncr->Size() != n) {
        if (stepIncr != 0) delete stepIncr;
        stepIncr = new Vector(n);
    }
    stepIncr->Zero();

    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
        opserr << "WarmStartNewton::solveCurrentStep() -";
        opserr << "the ConvergenceTest object failed in start()\n";
        return -3;
    }

    int result = -1;
    numIterations = 0;

    if (thePredictor != 0) {

        if (guessWork == 0 || guessWork->Size() != n) {
            if (guessWork != 0) delete guessWork;
            guessWork = new Vector(n);
        }
        Vector &dUguess = *guessWork;
        dUguess.Zero();

        if (thePredictor->predict(*theAnaModel, *theIntegrator, dUguess)) {

            bool accept = (thePredictor->confidence() >= minConfidence);

            if (accept) {
                double scale = this->referenceScale();
                if (scale > 0.0 && dUguess.Norm() > trustRadius * scale)
                    accept = false;
            }

            if (accept) {
                // First gate: never-worse against the residual we just
                // formed at the integrator's own predicted point.
                //
                // This is also what currently keeps an unconstrained
                // guess from corrupting ArcLength / DisplacementControl:
                // violating their constraint equation makes the residual
                // worse, so the guess is rejected. That is a guard, not
                // support for those integrators.
                double residBefore = theSOE->getB().Norm();

                if (theIntegrator->update(dUguess) < 0) {
                    opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
                    opserr << "the Integrator failed applying the predicted step\n";
                    return -4;
                }
                if (theIntegrator->formUnbalance() < 0) {
                    opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
                    opserr << "the Integrator failed in formUnbalance()\n";
                    return -2;
                }
                double residAfter = theSOE->getB().Norm();

                bool residualImproved = (residAfter < residBefore);
                bool lookaheadPassed = false;

                // Second gate: a single trial-point residual comparison
                // cannot tell whether the guess actually reduces the
                // number of iterations still needed -- it can look
                // locally better while leaving the state somewhere
                // Newton needs MORE total work to fix (overshoot, wrong-
                // direction curvature). Form the tangent and solve AT the
                // guessed point -- this is the mandatory first Newton
                // iteration regardless of accept/reject, not wasted work
                // -- and only commit to the guess if the resulting
                // correction is small relative to a typical full step.
                //
                // Compared against referenceScale() (the last converged
                // step's increment), NOT ||dUguess|| -- an extrapolation
                // guess is a SECOND difference of consecutive increments
                // and is intrinsically tiny relative to the increment
                // itself, so gating on the guess's own magnitude made the
                // bar nearly impossible to clear regardless of whether
                // the guess actually helped (empirically: 0% acceptance
                // up to factor=3, only 3.5% at factor=5, on a 20-story
                // frame where the ungated predictor got 42% acceptance
                // with a real ~8% iteration reduction on the accepted
                // steps). A leftover correction on the order of a full
                // step's typical size is the more sensible failure
                // signal: it means the guess bought us essentially
                // nothing, not that its own magnitude was outscaled.
                if (residualImproved) {
                    SOLUTION_ALGORITHM_tangentFlag = tangent;
                    if (theIntegrator->formTangent(tangent, iFactor, cFactor) < 0) {
                        opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
                        opserr << "the Integrator failed in formTangent()\n";
                        return -1;
                    }
                    if (theSOE->solve() < 0) {
                        opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
                        opserr << "the LinearSysOfEqn failed in solve()\n";
                        return -3;
                    }

                    double dXTrialNorm = theSOE->getX().Norm();
                    double scale = this->referenceScale();

                    lookaheadPassed = (scale <= 0.0) ||
                                       (dXTrialNorm <= lookaheadFactor * scale);

                    if (lookaheadPassed) {
                        // Commit: this solve IS the first Newton
                        // iteration, reused rather than repeated.
                        stepIncr->addVector(1.0, dUguess, 1.0);
                        stepIncr->addVector(1.0, theSOE->getX(), 1.0);

                        if (theIntegrator->update(theSOE->getX()) < 0) {
                            opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
                            opserr << "the Integrator failed in update()\n";
                            return -4;
                        }
                        if (theIntegrator->formUnbalance() < 0) {
                            opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
                            opserr << "the Integrator failed in formUnbalance()\n";
                            return -2;
                        }

                        result = theTest->test();
                        numIterations = 1;
                        this->record(numIterations);

                        usedPrediction = true;
                        numPredictedSteps++;
                    }
                }

                if (!usedPrediction) {
                    // Undo exactly, so the state entering the loop is
                    // identical to the no-predictor case. The lookahead
                    // solve (if it ran) does not need undoing -- it was
                    // never applied via update().
                    dUguess *= -1.0;
                    if (theIntegrator->update(dUguess) < 0) {
                        opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
                        opserr << "the Integrator failed reverting the predicted step\n";
                        return -4;
                    }
                    if (theIntegrator->formUnbalance() < 0) {
                        opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
                        opserr << "the Integrator failed in formUnbalance()\n";
                        return -2;
                    }
                    numRejectedGuesses++;
                }
            } else {
                numRejectedGuesses++;
            }
        }
    }

    if (!usedPrediction)
        numFallbackSteps++;

    // ------------------------------------------------------------------
    // NewtonRaphson iteration loop, unchanged, picking up from
    // numIterations (0 for a fallback/rejected step, 1 if the warm-start
    // guess's lookahead solve was just committed as iteration 1).
    // ------------------------------------------------------------------
    while (result == -1) {

        SOLUTION_ALGORITHM_tangentFlag = tangent;
        if (theIntegrator->formTangent(tangent, iFactor, cFactor) < 0) {
            opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
            opserr << "the Integrator failed in formTangent()\n";
            return -1;
        }

        if (theSOE->solve() < 0) {
            opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
            opserr << "the LinearSysOfEqn failed in solve()\n";
            return -3;
        }

        // Must be accumulated before the next solve() overwrites theSOE's
        // X vector. This is what closes the bootstrap loop: without it,
        // stepIncr only ever holds an accepted guess (never the Newton
        // loop's own corrections), so a step with no accepted guess -- the
        // only kind that can happen before any history exists -- reports a
        // ~zero converged increment to the predictor via observe(), which
        // then declines to ever predict (see ExtrapolationPredictor::
        // predict()'s n1/n2 near-zero guard). The predictor was being
        // starved of the exact signal it needs to start firing.
        stepIncr->addVector(1.0, theSOE->getX(), 1.0);

        if (theIntegrator->update(theSOE->getX()) < 0) {
            opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
            opserr << "the Integrator failed in update()\n";
            return -4;
        }

        if (theIntegrator->formUnbalance() < 0) {
            opserr << "WARNING WarmStartNewton::solveCurrentStep() -";
            opserr << "the Integrator failed in formUnbalance()\n";
            return -2;
        }

        result = theTest->test();
        numIterations++;
        this->record(numIterations);
    }

    if (result == -2) {
        opserr << "WarmStartNewton::solveCurrentStep() -";
        opserr << "the ConvergenceTest object failed in test()\n";
        return -3;
    }

    // ------------------------------------------------------------------
    // Bank diagnostics and feed the predictor.
    // ------------------------------------------------------------------
    if (usedPrediction)
        totalIterationsPredicted += numIterations;
    else
        totalIterationsFallback += numIterations;

    const Vector &convergedIncr = *stepIncr;

    if (lastConvergedIncr == 0 ||
        lastConvergedIncr->Size() != convergedIncr.Size()) {
        if (lastConvergedIncr != 0) delete lastConvergedIncr;
        lastConvergedIncr = new Vector(convergedIncr);
    } else {
        *lastConvergedIncr = convergedIncr;
    }

    if (thePredictor != 0)
        thePredictor->observe(*theAnaModel, *theIntegrator, convergedIncr);

    return result;
}

int
WarmStartNewton::domainChanged(void)
{
    if (stepIncr != 0) {
        delete stepIncr;
        stepIncr = 0;
    }
    if (lastConvergedIncr != 0) {
        delete lastConvergedIncr;
        lastConvergedIncr = 0;
    }
    if (guessWork != 0) {
        delete guessWork;
        guessWork = 0;
    }
    if (thePredictor != 0)
        thePredictor->reset();

    return 0;
}

int
WarmStartNewton::sendSelf(int cTag, Channel &theChannel)
{
    // An opaque NewtonPredictor cannot be serialized. Rather than
    // silently running Newton on the subdomains and something else on
    // the master, fail loudly. Parallel support needs a registration
    // tag and the predictor present on every process.
    if (thePredictor != 0) {
        opserr << "WarmStartNewton::sendSelf() - a predictor is attached; "
               << "parallel execution is not supported.\n";
        return -1;
    }

    static Vector data(6);
    data(0) = tangent;
    data(1) = iFactor;
    data(2) = cFactor;
    data(3) = trustRadius;
    data(4) = minConfidence;
    data(5) = lookaheadFactor;
    return theChannel.sendVector(this->getDbTag(), cTag, data);
}

int
WarmStartNewton::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static Vector data(6);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)
        return -1;

    tangent         = (int)data(0);
    iFactor         = data(1);
    cFactor         = data(2);
    trustRadius     = data(3);
    minConfidence   = data(4);
    lookaheadFactor = data(5);

    if (thePredictor != 0) {
        delete thePredictor;
        thePredictor = 0;
    }
    return 0;
}

void
WarmStartNewton::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) {
        s << "WarmStartNewton\n";
        s << "  predictor: "
          << (thePredictor != 0 ? thePredictor->getType()
                                : "none (equivalent to NewtonRaphson)") << "\n";
        s << "  trustRadius: " << trustRadius
          << "  minConfidence: " << minConfidence
          << "  lookaheadFactor: " << lookaheadFactor << "\n";
        s << "  steps warm-started: " << numPredictedSteps
          << "  fallback: " << numFallbackSteps
          << "  guesses rejected: " << numRejectedGuesses << "\n";
        s << "  avg iterations - warm-started: "
          << this->getAvgIterationsPredicted()
          << "  fallback: " << this->getAvgIterationsFallback() << "\n";
    }
}
