import 'dart:math';

import 'package:complex/complex.dart';

extension ComplexExt on Complex {

  double phase() {
    return atan2(imaginary, real);
  }

  double magnitude() {
    return sqrt(real * real + imaginary * imaginary);
  }
}

class ComplexVector {
  List<double> re;
  List<double> im;
  ComplexVector(this.re, this.im);
}
