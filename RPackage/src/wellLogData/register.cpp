#include <Rcpp.h>
#include <internal.h>
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
#include "basicFilterRPackage.h"
#include "fearnheadFilterRPackage.h"
#include "withoutReplacementRPackage.h"
extern "C" const char* package_name = "wellLogData";
R_CallMethodDef callMethods[] = 
{
	{"basicFilter", (DL_FUNC)&basicFilter, 11},
	{"fearnheadFilter", (DL_FUNC)&fearnheadFilter, 11},
	{"withoutReplacement", (DL_FUNC)&withoutReplacement, 11},
	{NULL, NULL, 0}
};
RcppExport void R_init_wellLogData(DllInfo *info)
{
	std::vector<R_CallMethodDef> callMethodsVector;
	R_CallMethodDef* mpMap2CallMethods = callMethods;
	while(mpMap2CallMethods->name != NULL) mpMap2CallMethods++;
	callMethodsVector.insert(callMethodsVector.begin(), callMethods, mpMap2CallMethods);

#ifdef CUSTOM_STATIC_RCPP
	R_CallMethodDef* RcppStartCallMethods = Rcpp_get_call();
	R_CallMethodDef* RcppEndCallMethods = RcppStartCallMethods;
	while(RcppEndCallMethods->name != NULL) RcppEndCallMethods++;
	callMethodsVector.insert(callMethodsVector.end(), RcppStartCallMethods, RcppEndCallMethods);
#endif
	R_CallMethodDef blank = {NULL, NULL, 0};
	callMethodsVector.push_back(blank);

	R_registerRoutines(info, NULL, &(callMethodsVector[0]), NULL, NULL);
#ifdef CUSTOM_STATIC_RCPP
	init_Rcpp_cache();
#endif
}
