/*
 *  DllHeader.h
 *
 *  Copyright 2011 Plexim GmbH. All rights reserved.
 *
 */

 #ifndef PLECSDLLHEADER_H_
 #define PLECSDLLHEADER_H_
 
 #if defined _WIN32
   #define DLLEXPORT __declspec(dllexport)
 #elif __GNUC__ >=4
   #define DLLEXPORT __attribute__ ((visibility("default")))
 #else
   #define DLLEXPORT
 #endif
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 #pragma pack(push, 4)
 struct SimulationSizes {
    int numInputs;      /**< The number of inputs that the DLL needs. */
    int numOutputs;     /**< The number of outputs that the DLL provides. */
    int numStates;      /**< The number of discrete states that the DLL needs. */
    int numParameters;  /**< The number of user parameters that the DLL needs. */
 };
 
 struct SimulationState {
    const double* const inputs;     /**< Array of input values (read-only). */
    double* const outputs;          /**< Array of output values (to fill by DLL). */
    double* const states;           /**< Array of discrete states (read/write). */
    const double* const parameters; /**< Array of parameters (read-only). */
    const double time;              /**< Current simulation time (read-only). */
    const char* errorMessage;       /**< Error message to set by DLL. */
    void* userData;                 /**< Pointer to any DLL data (untouched by PLECS). */
 };
 #pragma pack(pop)
 
 /** Required: DLL needs to set all fields in aSizes.
  * Called once before the simulation. */
 DLLEXPORT void plecsSetSizes(struct SimulationSizes* aSizes);
 
 /** Optional: DLL may acquire resources, initialize states and outputs.
  * Called once during the initialization of a new simulation. */
 DLLEXPORT void plecsStart(struct SimulationState* aState);
 
 /** Required: DLL needs to set outputs depending on inputs and states.
  * Called whenever the simulation time reaches a multiple of the sample Time. */
 DLLEXPORT void plecsOutput(struct SimulationState* aState);
 
 /** Optional: DLL may release any acquired resources.
  * Called when the simulation is finished, even when an error occured. */
 DLLEXPORT void plecsTerminate(struct SimulationState* aState);
 
 #ifdef __cplusplus
 }
 #endif
 
 #endif /* PLECSDLLHEADER_H_ */
 