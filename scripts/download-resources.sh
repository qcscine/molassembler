RDL_VERSION=1.1.2
OUTCOME_VERSION=2.1.4
JSON_VERSION=3.9.1

mkdir resources
wget -q -O resources/nauty27r1.tar.gz http://pallini.di.uniroma1.it/nauty27r1.tar.gz
wget -q -O resources/v${RDL_VERSION}.tar.gz https://github.com/rareylab/RingDecomposerLib/archive/v${RDL_VERSION}.tar.gz
mkdir resources/outcome
wget -q -O resources/outcome/outcome.hpp https://raw.githubusercontent.com/ned14/outcome/v${OUTCOME_VERSION}/single-header/outcome.hpp
wget -q -O resources/outcome/LICENSE https://raw.githubusercontent.com/ned14/outcome/v${OUTCOME_VERSION}/Licence.txt
mkdir resources/nlohmann
wget -q -O resources/nlohmann/json.hpp https://raw.githubusercontent.com/nlohmann/json/v${JSON_VERSION}/single_include/nlohmann/json.hpp
wget -q -O resources/nlohmann/LICENSE https://raw.githubusercontent.com/nlohmann/json/v${JSON_VERSION}/LICENSE.MIT 
