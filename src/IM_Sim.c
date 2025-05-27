/// \file		psim.c
///
/// \brief		This dll is unique to PSIM and Plecs. 
///
/// \author		Uriel Abe Contardi (urielcontardi@hotmail.com)
/// \date		27-02-2025
///
/// \version	1.0
///
/// \note		Revisions:
/// 			27-02-2025 <urielcontardi@hotmail.com>
/// 			First revision.
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               INCLUDES                                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include <DllHeader.h>
#include <windows.h>
#include <IM_Model.h>

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           DEFINES AND MACROS                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                      LOCAL TYPEDEFS AND STRUCTURES                       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
typedef enum
{
	VA=0,
    VB,
    VC,
    TLOAD,
    RS,
    RR,
    LM,
    LS,
    LR,
    J,
    NPP,
    TS,
    MODEL,
	TOTAL_INPUTS
} Inputs_t;

typedef enum
{
    IA=0,
    IB,
    IC,
    WR,
    WMEC,
    TE,
    VALPHA,
    VBETA,
	TOTAL_OUTPUTS
} Outputs_t;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                        LOCAL FUNCTIONS PROTOTYPES                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
static IM_Model_t* _getIMInstance(void);
static void _setParameters(IM_Model_t* model, const double *simInputs);
static void _setInputs(IM_Model_t* model, const double *simInputs);
static void _runSimulation(IM_Model_t* pIM, double t, double delt, const double* simInputs, double* simOut);

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                      STATIC VARIABLES AND CONSTANTS                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            EXPORTED FUNCTIONS                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// PSIM Calback function
__declspec(dllexport) void simuser(double t, double delt, double *simInputs, double *simOut);

__declspec(dllexport) void simuser(double t, double delt, double* simInputs, double* simOut)
{
    IM_Model_t* pIM = _getIMInstance();
    _runSimulation(pIM, t, delt, simInputs, simOut);
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            EXPORTED FUNCTIONS                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// PLECS Calback functions
DLLEXPORT void plecsSetSizes(struct SimulationSizes* aSizes)
{
   aSizes->numInputs = TOTAL_INPUTS;
   aSizes->numOutputs = TOTAL_OUTPUTS;
   aSizes->numStates = 0;
   aSizes->numParameters = 0;
}

//This function is automatically called at the beginning of the simulation
DLLEXPORT void plecsStart(struct SimulationState* aState)
{

}

//This function is automatically called every sample time
//output is written to DLL output port after the output delay
DLLEXPORT void plecsOutput(struct SimulationState* aState)
{
    const double* simInputs = aState->inputs;
    double* simOut = aState->outputs;
    double t = aState->time;
    const double DELT = 1e-7;

    IM_Model_t* pIM = _getIMInstance();
    _runSimulation(pIM, t, DELT, simInputs, simOut);
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              LOCAL FUNCTIONS                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
static void _runSimulation(IM_Model_t* pIM, double t, double delt, const double* simInputs, double* simOut)
{
    if (t <= delt) {
        IM_Init(pIM);
        _setParameters(pIM, simInputs);
        IM_TypeModel(pIM, (IMType)simInputs[MODEL]);
    }

    IM_SimulateStep(pIM);
    _setInputs(pIM, simInputs);

    // Copy outputs to simOut
    simOut[IA]   = (double)pIM->out.ia;
    simOut[IB]   = (double)pIM->out.ib;
    simOut[IC]   = (double)pIM->out.ic;
    simOut[WR]   = (double)pIM->out.wr;
    simOut[WMEC] = (double)pIM->out.wmec;
    simOut[TE]   = (double)pIM->out.Te;

    double* privData = (double*)pIM->priv;
    simOut[VALPHA] = privData[0];
    simOut[VBETA]  = privData[1];
}

static IM_Model_t* _getIMInstance(void)
{
	static IM_Model_t model;
	return &model;
}

static void _setParameters(IM_Model_t* model, const double *simInputs)
{
    IMParams params;
    params.Rs = simInputs[RS];
    params.Rr = simInputs[RR];
    params.Lm = simInputs[LM];
    params.Ls = simInputs[LS];
    params.Lr = simInputs[LR];
    params.J = simInputs[J];
    params.npp = simInputs[NPP];
    params.Ts = simInputs[TS];
    IM_SetParams(model, &params);
}

static void _setInputs(IM_Model_t* model,const double *simInputs)
{
    IMInputs inp;
    inp.Va = simInputs[VA];
    inp.Vb = simInputs[VB];
    inp.Vc = simInputs[VC];
    inp.Tload = simInputs[TLOAD];
    IM_SetInputs(model, &inp);
}
