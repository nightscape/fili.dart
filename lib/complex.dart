import 'dart:math';

class Complex {
  double re;
  double im;
  Complex(this.re, this.im);
  Complex div(Complex q) {
    var a = re;
    var b = im;
    var c = q.re;
    var d = q.im;
    var n = c * c + d * d;
    return new Complex(
      (a * c + b * d) / n,
      (b * c - a * d) / n,
    );
  }

  Complex mul(Complex q) {
    var a = re;
    var b = im;
    var c = q.re;
    var d = q.im;
    return new Complex(
      a * c - b * d,
      (a + b) * (c + d) - a * c - b * d,
    );
  }

  Complex add(Complex q) {
    return new Complex(
      re + q.re,
      im + q.im,
    );
  }

  Complex sub(Complex q) {
    return new Complex(
      re - q.re,
      im - q.im,
    );
  }

  double phase() {
    return atan2(im, re);
  }

  double magnitude() {
    return sqrt(re * re + im * im);
  }
}

class ComplexVector {
  List<double> re;
  List<double> im;
  ComplexVector(this.re, this.im);
}
