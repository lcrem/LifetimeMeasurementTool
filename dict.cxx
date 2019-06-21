// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "UsefulFunctions.h"
#include "RawDigitiser.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_RawDigitiser(void *p = 0);
   static void *newArray_RawDigitiser(Long_t size, void *p);
   static void delete_RawDigitiser(void *p);
   static void deleteArray_RawDigitiser(void *p);
   static void destruct_RawDigitiser(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RawDigitiser*)
   {
      ::RawDigitiser *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RawDigitiser >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RawDigitiser", ::RawDigitiser::Class_Version(), "RawDigitiser.h", 12,
                  typeid(::RawDigitiser), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RawDigitiser::Dictionary, isa_proxy, 4,
                  sizeof(::RawDigitiser) );
      instance.SetNew(&new_RawDigitiser);
      instance.SetNewArray(&newArray_RawDigitiser);
      instance.SetDelete(&delete_RawDigitiser);
      instance.SetDeleteArray(&deleteArray_RawDigitiser);
      instance.SetDestructor(&destruct_RawDigitiser);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RawDigitiser*)
   {
      return GenerateInitInstanceLocal((::RawDigitiser*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RawDigitiser*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RawDigitiser::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RawDigitiser::Class_Name()
{
   return "RawDigitiser";
}

//______________________________________________________________________________
const char *RawDigitiser::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RawDigitiser*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RawDigitiser::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RawDigitiser*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RawDigitiser::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RawDigitiser*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RawDigitiser::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RawDigitiser*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RawDigitiser::Streamer(TBuffer &R__b)
{
   // Stream an object of class RawDigitiser.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RawDigitiser::Class(),this);
   } else {
      R__b.WriteClassBuffer(RawDigitiser::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RawDigitiser(void *p) {
      return  p ? new(p) ::RawDigitiser : new ::RawDigitiser;
   }
   static void *newArray_RawDigitiser(Long_t nElements, void *p) {
      return p ? new(p) ::RawDigitiser[nElements] : new ::RawDigitiser[nElements];
   }
   // Wrapper around operator delete
   static void delete_RawDigitiser(void *p) {
      delete ((::RawDigitiser*)p);
   }
   static void deleteArray_RawDigitiser(void *p) {
      delete [] ((::RawDigitiser*)p);
   }
   static void destruct_RawDigitiser(void *p) {
      typedef ::RawDigitiser current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RawDigitiser

namespace {
  void TriggerDictionaryInitialization_dict_Impl() {
    static const char* headers[] = {
"UsefulFunctions.h",
"RawDigitiser.h",
0
    };
    static const char* includePaths[] = {
"/home/lindac/DUNE/rootInstall/root/include",
"/home/lindac/DUNE/utils//include",
"/home/lindac/DUNE/rootInstall/root/include",
"/home/lindac/DUNE/LifetimeMeasurementTool/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$RawDigitiser.h")))  RawDigitiser;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "UsefulFunctions.h"
#include "RawDigitiser.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RawDigitiser", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dict() {
  TriggerDictionaryInitialization_dict_Impl();
}
