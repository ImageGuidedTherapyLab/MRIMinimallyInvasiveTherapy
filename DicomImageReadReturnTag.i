%module DicomImageReadReturnTag
%include "std_string.i"
%{
/* Put headers and other declarations here */
std::string GetDicomTag(const std::string &InputFile, const std::string &TagID, const std::string &DicomDictionary);
%}

std::string GetDicomTag(const std::string &InputFile, const std::string &TagID, const std::string &DicomDictionary);

