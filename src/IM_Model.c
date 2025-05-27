/// \file		IM_Model.c
///
/// \brief	    
///
/// \author		Uriel Abe Contardi (urielcontardi@hotmail.com)
/// \date		26-01-2025
///
/// \version	1.0
///
/// \note		Revisions:
/// 			26-01-2025 <urielcontardi@hotmail.com>
/// 			First revision.
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               INCLUDES                                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include "IM_Model.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           DEFINES AND MACROS                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#define M_PI 3.14159265358979323846

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                      LOCAL TYPEDEFS AND STRUCTURES                       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
typedef struct {
    double is_alpha;
    double is_beta;
    double ir_alpha;
    double ir_beta;
    double fluxR_alpha;
    double fluxR_beta;
    double wr;
    double wm;
    double Te;
    double isd;
    double isq;
    double fluxRd;
    double angleR;
} IM_States_t;

typedef struct {
    double valpha;
    double vbeta;
    double v0;
} IM_InternalInputs_t;

typedef struct {
    IM_InternalInputs_t inp;
    IM_States_t out;
} IM_PrivateData_t;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                        LOCAL FUNCTIONS PROTOTYPES                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
static void _vabc2AphaBeta(IM_Model_t *self);
static void _updateOutputs(IM_Model_t *self);

static void _updateStatesA(IM_Model_t *self);
static void _updateStatesB(IM_Model_t *self);
static void _updateStatesB2(IM_Model_t *self);
static void _updateStatesC(IM_Model_t *self);
static void _updateStatesD(IM_Model_t *self);

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                      STATIC VARIABLES AND CONSTANTS                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

static IM_PrivateData_t _privateData;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            EXPORTED FUNCTIONS                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void IM_Init(IM_Model_t *self) {
    if (self == NULL) {
        return; 
    }

    // Clear the structure
    memset(self, 0, sizeof(IM_Model_t));

    // Allocation of the for the priv data
    self->priv = (void*)&_privateData;
    if (self->priv == NULL) {
        return;
    }

    // Clear priv data
    memset(self->priv, 0, sizeof(IM_PrivateData_t));
}

void IM_SetParams(IM_Model_t *self, const IMParams *params) {
    if (self == NULL || params == NULL) {
        return;
    }

    memcpy(&self->params, params, sizeof(IMParams));
}

void IM_SetInputs(IM_Model_t *self, const IMInputs *inputs) {
    if (self == NULL || inputs == NULL) {
        return;
    }

    memcpy(&self->inp, inputs, sizeof(IMInputs));
    _vabc2AphaBeta(self);
}

void IM_TypeModel(IM_Model_t *self, IMType model)
{
    self->type = model;
    memset(self->priv, 0, sizeof(IM_PrivateData_t));
}

void IM_SimulateStep(IM_Model_t *self) {

    // Compute selected model
    switch (self->type )
    {
        default:
        case MODEL_A:
            _updateStatesA(self);
            break;

        case MODEL_B:
            _updateStatesB(self);
            break;

        case MODEL_B2:
            _updateStatesB2(self);
            break;

        case MODEL_C:
            _updateStatesC(self);
            break;

        case MODEL_D:
            _updateStatesD(self);
            break;
    
    }

    _updateOutputs(self);
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              LOCAL FUNCTIONS                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void _vabc2AphaBeta(IM_Model_t *self) {

    // Get the input values
    double Va = self->inp.Va;
    double Vb = self->inp.Vb;
    double Vc = self->inp.Vc;
    
    // Clark Transform (abc -> αβ)
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_InternalInputs_t *intInputs = &privateData->inp;
    intInputs->valpha = (2.0 / 3.0) * (Va - 0.5 * Vb - 0.5 * Vc);
    intInputs->vbeta =  (1.0 / sqrt(3.0)) * (Vb - Vc);
    intInputs->v0 = (1.0/3.0) * (Va+Vb+Vc);

}

void _updateStatesA(IM_Model_t *self) {
    /// Article: FPGA based real-time model of three-phase induction machine
    /// https://ieeexplore.ieee.org/document/8395534
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // MIT Parameters
    double Rs  = self->params.Rs;
    double Rr  = self->params.Rr;
    double Lm  = self->params.Lm;
    double Ls  = self->params.Ls + self->params.Lm;
    double Lr  = self->params.Lr + self->params.Lm;
    double J   = self->params.J; 
    double npp = self->params.npp;
    double Ts  = self->params.Ts;

    // Inputs (v_alpha, v_beta)
    double v_alpha_s = intInputs->valpha;
    double v_beta_s  = intInputs->vbeta;
    double Tload     = self->inp.Tload;

    // Article Equations :
    double dfluxR_alpha = (Rr*Lm/Lr)*out->is_alpha - out->wr * out->fluxR_beta - (Rr/Lr) * out->fluxR_alpha;
    double dfluxR_beta = (Rr*Lm/Lr)*out->is_beta + out->wr * out->fluxR_alpha - (Rr/Lr) * out->fluxR_beta;
    double dis_alpha = (Lr/(Lm*Lm - Lr*Ls))*(Rs * out->is_alpha + (Lm/Lr)*dfluxR_alpha - v_alpha_s);
    double dis_beta = (Lr/(Lm*Lm - Lr*Ls))*(Rs * out->is_beta + (Lm/Lr)*dfluxR_beta - v_beta_s);

    double Te =  (3.0f/2.0f)*(npp*Lm/Lr)*(out->fluxR_alpha*out->is_beta - out->fluxR_beta*out->is_alpha);
    double dwm = (Te - Tload)/J;

    // Compute next value
    out->is_alpha    = out->is_alpha + dis_alpha * Ts;
    out->is_beta     = out->is_beta + dis_beta * Ts;
    out->fluxR_alpha = out->fluxR_alpha + dfluxR_alpha * Ts;
    out->fluxR_beta  = out->fluxR_beta + dfluxR_beta * Ts;
    out->wm          = out->wm + dwm * Ts;
    out->wr          = out->wm * npp;
    out->Te          = Te;

}

void _updateStatesB(IM_Model_t *self) {
    /// Article: Real-Time Emulator of an Induction Motor: FPGA-based Implementation
    /// https://ieeexplore.ieee.org/document/6421152/
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // MIT Parameters
    double Rs  = self->params.Rs;
    double Rr  = self->params.Rr;
    double Lm  = self->params.Lm;
    double Ls  = self->params.Ls + self->params.Lm;
    double Lr  = self->params.Lr + self->params.Lm;
    double J   = self->params.J; 
    double npp = self->params.npp;
    double Ts  = self->params.Ts;

    // Inputs (v_alpha, v_beta)
    double v_alpha_s = intInputs->valpha;
    double v_beta_s  = intInputs->vbeta;
    double Tload     = self->inp.Tload;

    // Constanst
    double Tr = Lr/Rr;
    double sigma = 1.0 - ((Lm*Lm)/(Ls*Lr));
    double K = Lm/(sigma * Ls * Lr);
    double gama = (Rs/(sigma * Ls)) + ((Rr*Lm*Lm)/(sigma*Ls*Lr*Lr));
    
    // Derivative
    double dis_alpha = -gama*out->is_alpha + (out->fluxR_alpha*K/Tr) + npp*out->wm*K*out->fluxR_beta + (v_alpha_s/(sigma*Ls));
    double dis_beta = -gama*out->is_beta + (out->fluxR_beta*K/Tr) - npp*out->wm*K*out->fluxR_alpha + (v_beta_s/(sigma*Ls));
    double dfluxR_alpha = (out->is_alpha*Lm/Tr) - out->fluxR_alpha/Tr - npp*out->wm*out->fluxR_beta;
    double dfluxR_beta = (out->is_beta*Lm/Tr) - out->fluxR_beta/Tr + npp*out->wm*out->fluxR_alpha;

    // Note: The torque equation is different from the article; comparing with other articles, I noticed the absence of a 3/2 factor.
    // By adding the 3/2 factor, the model behaves much closer to the PSIM model. I will check this further later.
    double dwm = (3.0f/2.0f)*(npp*Lm/(J*Lr))*(out->fluxR_alpha*out->is_beta - out->fluxR_beta*out->is_alpha) - (Tload/J);

    // Compute next value
    out->is_alpha = out->is_alpha + dis_alpha * Ts;
    out->is_beta  = out->is_beta + dis_beta * Ts;
    out->fluxR_alpha = out->fluxR_alpha + dfluxR_alpha * Ts;
    out->fluxR_beta  = out->fluxR_beta + dfluxR_beta * Ts;
    out->wm = out->wm + dwm * Ts;

}

void _updateStatesB2(IM_Model_t *self) {
    /// Based on Article: Real-Time Emulator of an Induction Motor: FPGA-based Implementation
    /// https://ieeexplore.ieee.org/document/6421152/
    // Here, the same algorithm as _updateStatesB was implemented, but without using the auxiliary constants
    // used by the author of the article.
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // MIT Parameters
    double Rs  = self->params.Rs;
    double Rr  = self->params.Rr;
    double Lm  = self->params.Lm;
    double Ls  = self->params.Ls + self->params.Lm;
    double Lr  = self->params.Lr + self->params.Lm;
    double J   = self->params.J; 
    double npp = self->params.npp;
    double Ts  = self->params.Ts;

    // Inputs (v_alpha, v_beta)
    double v_alpha_s = intInputs->valpha;
    double v_beta_s  = intInputs->vbeta;
    double Tload     = self->inp.Tload;

    // Constanst
    double k = 1.0f / (Lr * (Lm * Lm - Lr * Ls));
    
    // Derivative
    double dis_alpha = k * (Lm*Lm*Rr*out->is_alpha - Lm*Lr*npp*out->wm*out->fluxR_beta - Lm*Rr*out->fluxR_alpha + Lr*Lr*Rs*out->is_alpha - Lr*Lr*v_alpha_s);
    double dis_beta = k * (Lm*Lm*Rr*out->is_beta + Lm*Lr*npp*out->wm*out->fluxR_alpha - Lm*Rr*out->fluxR_beta + Lr*Lr*Rs*out->is_beta - Lr*Lr*v_beta_s);
    double dfluxR_alpha = (Lm*Rr*out->is_alpha - Lr*npp*out->wm*out->fluxR_beta - Rr*out->fluxR_alpha)/Lr;
    double dfluxR_beta = (Lm*Rr*out->is_beta + Lr*npp*out->wm*out->fluxR_alpha - Rr*out->fluxR_beta)/Lr;

    double Te = (3.0f/2.0f)*(npp*Lm/Lr)*(out->fluxR_alpha*out->is_beta - out->fluxR_beta*out->is_alpha);
    double dwm = (Te - Tload)/J;

    // Compute next value
    out->is_alpha    = out->is_alpha + dis_alpha * Ts;
    out->is_beta     = out->is_beta + dis_beta * Ts;
    out->fluxR_alpha = out->fluxR_alpha + dfluxR_alpha * Ts;
    out->fluxR_beta  = out->fluxR_beta + dfluxR_beta * Ts;
    out->wm          = out->wm + dwm * Ts;
    out->wr          = out->wm * npp;
    out->Te          = Te;

}

void _updateStatesC(IM_Model_t *self) {
    /// Article: C++ based dynamic model of AC induction motor in discrete time domain
    /// https://ieeexplore.ieee.org/document/7512961
    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    // MIT Parameters
    double Rs  = self->params.Rs;
    double Rr  = self->params.Rr;
    double Lm  = self->params.Lm;
    double Ls  = self->params.Ls + self->params.Lm;
    double Lr  = self->params.Lr + self->params.Lm;
    double J   = self->params.J; 
    double npp = self->params.npp;
    double Ts  = self->params.Ts;

    // Inputs (v_alpha, v_beta)
    double v_alpha_s = intInputs->valpha;
    double v_beta_s  = intInputs->vbeta;
    double Tload     = self->inp.Tload;

    // Constanst
    double k = 1.0f / (Lm * Lm - Lr * Ls);
    
    // Derivative
    double dis_alpha = k * (-Lm*Lm*out->wr*out->is_beta - Lm*Lr*out->wr*out->ir_beta - Lm*Rr*out->ir_alpha + Lr*Rs*out->is_alpha - Lr*v_alpha_s);
    double dis_beta  = k * (Lm*Lm*out->wr*out->is_alpha + Lm*Lr*out->wr*out->ir_alpha - Lm*Rr*out->ir_beta + Lr*Rs*out->is_beta - Lr*v_beta_s);
    double dir_alpha = k * (Lm*Ls*out->wr*out->is_beta - Lm*Rs*out->is_alpha + Lm*v_alpha_s + Lr*Ls*out->wr*out->ir_beta + Ls*Rr*out->ir_alpha);
    double dir_beta  = k * (-Lm*Ls*out->wr*out->is_alpha - Lm*Rs*out->is_beta + Lm*v_beta_s - Lr*Ls*out->wr*out->ir_alpha + Ls*Rr*out->ir_beta);

    // Update Torque / Speed
    double Te = (3.0f/2.0f)* npp * Lm * (out->is_beta * out->ir_alpha - out->is_alpha * out->ir_beta);
    double dwm = (Te - Tload)/J;

    // Update Output
    out->is_alpha = out->is_alpha + dis_alpha * Ts;
    out->is_beta  = out->is_beta  + dis_beta * Ts;
    out->ir_alpha = out->ir_alpha + dir_alpha * Ts;
    out->ir_beta  = out->ir_beta  + dir_beta * Ts;
    out->wm       = out->wm + Ts * dwm;
    out->wr       = out->wm * npp;
    out->Te       = Te;
}

void _updateStatesD(IM_Model_t *self) {

    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
//    IM_InternalInputs_t *intInputs = &privateData->inp;

    // MIT Parameters
    double Rs  = self->params.Rs;
    double Rr  = self->params.Rr;
    double Lm  = self->params.Lm;
    double Ls  = self->params.Ls + self->params.Lm;
    double Lr  = self->params.Lr + self->params.Lm;
    double J   = self->params.J; 
    double npp = self->params.npp;
    double Ts  = self->params.Ts;

    // Inverse Gamma parameters
    double gamma = Lm/(Lr+Lm);
    double Lsigma = Ls + gamma * Lr;
    double LM = gamma * Lm;
    double RR = gamma * gamma * Rr;
    double Rsigma = Rs + RR;
    double a = RR/LM;

    // Inputs
    double Tload     = self->inp.Tload;
    
    // abc -> dq
    double theta = out->wr;
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double cos_theta_120 = cos(theta - M_PI * 2 / 3); // cos(θ - 120°)
    double sin_theta_120 = sin(theta - M_PI * 2 / 3); // sin(θ - 120°)
    double cos_theta_240 = cos(theta + M_PI * 2 / 3); // cos(θ + 120°)
    double sin_theta_240 = sin(theta + M_PI * 2 / 3); // sin(θ + 120°)

    double Va = self->inp.Va;
    double Vb = self->inp.Vb;
    double Vc = self->inp.Vc;
    
    float Vsd = (2.0 / 3.0) * (Va * cos_theta + Vb * cos_theta_120 + Vc * cos_theta_240);
    float Vsq = (2.0 / 3.0) * (-Va * sin_theta - Vb * sin_theta_120 - Vc * sin_theta_240);
    
    // Currents
    float d_isd = 0.0f;
    float d_isq = 0.0f;
    d_isd = (Ts/Lsigma)*(Vsd - Rsigma*out->isd + Lsigma*out->wr*out->isq + a*out->fluxRd);
    d_isq = (Ts/Lsigma)*(Vsq - Rsigma*out->isq - Lsigma*out->wr*out->isd - npp*out->wm*out->fluxRd);
    
    // Flux
    float d_fluxRd = 0.0f;
    float d_angleR = 0.0f;
    d_fluxRd = Ts * (RR * out->isd + a * out->fluxRd);
    d_angleR = Ts * (RR * out->isq/out->fluxRd) + npp * out->wm;
    
    // Torque
    out->Te = (3.0/2.0) * npp * out->fluxRd * out->isq;
    
    // Wm
    float d_wm = 0.0f;
    d_wm = Ts * (out->Te - Tload) / J;

    // Update states
    out->wr = (RR * out->isq / out->fluxRd) + npp * out->wm;
    out->wm += d_wm;
    out->isd += d_isd;
    out->isq += d_isq;
    out->fluxRd += d_fluxRd;
    out->angleR += d_angleR;

}

void _updateOutputs(IM_Model_t *self) {

    IM_PrivateData_t *privateData = (IM_PrivateData_t *)self->priv;
    IM_States_t *out = &privateData->out;
    IM_InternalInputs_t *intInputs = &privateData->inp;

    if (self->type == MODEL_D)
    {
        double theta = out->wr;
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);
        double cos_theta_120 = cos(theta - M_PI * 2 / 3); // cos(θ - 120°)
        double sin_theta_120 = sin(theta - M_PI * 2 / 3); // sin(θ - 120°)
        double cos_theta_240 = cos(theta + M_PI * 2 / 3); // cos(θ + 120°)
        double sin_theta_240 = sin(theta + M_PI * 2 / 3); // sin(θ + 120°)
    
        // dq currents
        double isd = out->isd;
        double isq = out->isq;
    
        // dq -> abc
        self->out.ia = cos_theta * isd - sin_theta * isq;
        self->out.ib = cos_theta_120 * isd - sin_theta_120 * isq;
        self->out.ic = cos_theta_240 * isd - sin_theta_240 * isq;
    } else {
        self->out.ia = out->is_alpha + intInputs->v0;
        self->out.ib = -0.5 * out->is_alpha + (sqrt(3.0) / 2.0) * out->is_beta + intInputs->v0;
        self->out.ic = -0.5 * out->is_alpha - (sqrt(3.0) / 2.0) * out->is_beta + intInputs->v0;
    }

    self->out.wmec = out->wm;
    self->out.wr = out->wr;
    self->out.Te = out->Te;

}
