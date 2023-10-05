#include <HE_word.hpp>

helib::SecKey * HE_word::sk = NULL;
bool HE_word::verfier_mode = false;

HE_word::HE_word(/* args */)
{
}

HE_word::HE_word(HE_word &other){
  elemCT = other.elemCT;
  elemHE = new helib::Ctxt(*(other.elemHE));
  cleartext = other.cleartext;
}

HE_word::HE_word(uint64_t x){
  elemCT = x;
}

HE_word::HE_word(helib::Ctxt * x)
{
  elemHE = x;
  cleartext = false;
}

HE_word::~HE_word(){
  delete elemHE;
}
