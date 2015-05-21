%module itkUtilities
%include "std_string.i"
%include "exception.i"
/*
  Possible Swig Exceptions

     http://www.swig.org/Doc1.3/Library.html#Library_nn16

            SWIG_MemoryError
            SWIG_IOError
            SWIG_RuntimeError
            SWIG_IndexError
            SWIG_TypeError
            SWIG_DivisionByZero
            SWIG_OverflowError
            SWIG_SyntaxError
            SWIG_ValueError
            SWIG_SystemError

  Possible std::exception

     http://www.aoc.nrao.edu/php/tjuerges/ALMA/STL/html-3.4.6/classstd_1_1exception.html

      * std::exception
         * std::bad_alloc
         * std::bad_cast
         * std::bad_exception
         * std::bad_typeid
         * std::ios_base::failure
         * std::logic_error
               o std::domain_error
               o std::invalid_argument
               o std::length_error
               o std::out_of_range 
         * std::runtime_error
               o std::overflow_error
               o std::range_error
               o std::underflow_error 
*/
%exception {
  try {
    $action
  } catch (const itk::ExceptionObject& e) {
    SWIG_exception(SWIG_IOError, e.what());
  } catch (const std::logic_error& e) {
    SWIG_exception(SWIG_ValueError, e.what());
  } catch (...) {
    SWIG_exception(SWIG_RuntimeError, "Don't Know What Happened...");
  }
}
%include "std_vector.i"
/* Instantiate templates used by example */
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}

%{
/* Put headers and other declarations here */
#include "itkImage.h"
std::string GetDicomTag(const std::string &InputFile, const std::string &TagID, const std::string &DicomDictionary);
std::vector<double>  GetPixelValue(const std::string &InputFile, 
                                   const std::vector<int> inputIndex );
%}

std::string GetDicomTag(const std::string &InputFile, const std::string &TagID, const std::string &DicomDictionary);

std::vector<double>  GetPixelValue(const std::string &InputFile, 
                                   const std::vector<int> inputIndex );
