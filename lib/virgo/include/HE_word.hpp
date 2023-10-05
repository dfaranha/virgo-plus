#pragma once
#include <helib/helib.h>
#include <iostream>

class HE_word
{
private:
  helib::Ctxt * elemHE;
  uint64_t elemCT;
  bool cleartext = true;
public:
  static bool verfier_mode;
  static helib::SecKey * sk;
  HE_word(/* args */);
  HE_word(HE_word &other);
  HE_word(uint64_t x);
  HE_word(helib::Ctxt * x);
  ~HE_word();

  HE_word operator+(const HE_word &other) const;
  HE_word operator-(const HE_word &other) const;
  HE_word operator*(const HE_word &other) const;
  HE_word operator%(const HE_word &other) const;
};

