import 'dart:math';

extension MathExt on double {
  double sinh() => (exp(this) - exp(-this)) / 2;
  double sinc() => sin(pi * this) / (pi * this);
}
