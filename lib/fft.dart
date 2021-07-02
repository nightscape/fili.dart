import 'dart:math';

import 'package:fili.dart/complex.dart';
import 'package:fili.dart/math_ext.dart';

class Foo {
  double Function(double n, double N, double a) calc;
  List<double> values;
  double correction;
  Foo({required this.calc, required this.values, required this.correction});
}

class Fft {
  int radix;
  late FftData fft;
  late Map<String, Foo> windowCalculation;
  /*
    windowCalculation: {
        [String key]: { calc: (...List<dynamic> params) => number, List<double> values;, double correction; }
        rectangular: { calc: () => number; values: never[]; double correction;; }; none: { calc: () => number; values: never[]; double correction;; }; hanning: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hamming: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; tukery: { calc: (dynamic double n, double N, double a, a) => number; values: never[]; double correction;; }; cosine: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; lanczos: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; triangular: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; bartlett: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; gaussian: { calc: (dynamic double n, double N, double a, a) => number; values: never[]; double correction;; }; bartlettHanning: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; blackman: { calc: (dynamic double n, double N, double a, a) => number; values: never[]; double correction;; }; blackmanHarris: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; nuttall3: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; nuttall3a: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; nuttall3b: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; nuttall4: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; nuttall4a: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; nuttall4b: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; nuttall4c: { calc: (double n, double N, double a) => number; values: never[]; double correction;; };
        // fast decaying flat top
        sft3f: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; sft4f: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; sft5f: { calc: (double n, double N, double a) => number; values: never[]; double correction;; };
        // minimum sidelobe flat top
        sft3m: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; sft4m: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; sft5m: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; nift: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hpft: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; srft: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hft70: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hft95: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hft90d: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hft116d: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hft144d: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hft196d: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hft223d: { calc: (double n, double N, double a) => number; values: never[]; double correction;; }; hft248d: { calc: (double n, double N, double a) => number; values: never[]; double correction;; };
    };
    */
  static bool isPowerOfTwo(int value) {
    final one = value & (value - 1);
    final two = ~one;

    return (value & value - 1) == 0;
  }

  Fft(this.radix) {
    if (!Fft.isPowerOfTwo(radix)) {
      throw new Exception('Radix can only be a power of 2');
    }
    this.fft = FftData(
      length: radix,
      buffer: List.filled(radix, 0.0),
      re: List.filled(radix, 0.0),
      im: List.filled(radix, 0.0),
      reI: List.filled(radix, 0.0),
      imI: List.filled(radix, 0.0),
      twiddle: List.filled(radix, 0),
      sinTable: List.filled(radix - 1, 0.0),
      cosTable: List.filled(radix - 1, 0.0),
    );
    var TPI = 2 * pi;
    var bits = (log(radix) / ln2).floor();
    for (var i = this.fft.sinTable.length - 1; i >= 0; i--) {
      this.fft.sinTable[i] = sin(TPI * (i / radix));
      this.fft.cosTable[i] = cos(TPI * (i / radix));
    }
    var nh = radix >> 1;
    var i = 0;
    var j = 0;
    for (;;) {
      this.fft.twiddle[i] = j;
      if (++i >= radix) {
        break;
      }
      bits = nh;
      while (bits <= j) {
        j -= bits;
        bits >>= 1;
      }
      j += bits;
    }
    // good explanation in https://holometer.fnal.gov/GH_FFT.pdf
    var PI = pi;
    var PI2 = pi * 2;
    var E = e;
    this.windowCalculation = {
      "rectangular": Foo(
        calc: (double n, double N, double a) {
          return 1;
        },
        values: [],
        correction: 1,
      ),
      "none": Foo(
        calc: (double n, double N, double a) {
          return 1;
        },
        values: [],
        correction: 1,
      ),
      "hanning": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.5 * (1 - cos(z));
        },
        values: [],
        correction: 2,
      ),
      "hamming": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.54 - 0.46 * cos(z);
        },
        values: [],
        correction: 1.8518999946875638,
      ),
      "tukery": Foo(
        calc: (double n, double N, double a) {
          if (n < (a * (N - 1)) / 2) {
            return 0.5 * (1 + cos(PI * (((2 * n) / (a * (N - 1))) - 1)));
          } else if ((N - 1) * (1 - (a / 2)) < n) {
            return 0.5 *
                (1 + cos(PI * (((2 * n) / (a * (N - 1))) - (2 / a) + 1)));
          } else {
            return 1;
          }
        },
        values: [],
        correction: 4 / 3,
      ),
      "cosine": Foo(
        calc: (double n, double N, double a) {
          return sin((PI * n) / (N - 1));
        },
        values: [],
        correction: 1.570844266360796,
      ),
      "lanczos": Foo(
        calc: (double n, double N, double a) {
          return (((2 * n) / (N - 1)) - 1).sinc();
        },
        values: [],
        correction: 1.6964337576195783,
      ),
      "triangular": Foo(
        calc: (double n, double N, double a) {
          return (2 / (N + 1)) * (((N + 1) / 2) - (n - ((N - 1) / 2)).abs());
        },
        values: [],
        correction: 2,
      ),
      "bartlett": Foo(
        calc: (double n, double N, double a) {
          return (2 / (N - 1)) * (((N - 1) / 2) - (n - ((N - 1) / 2)).abs());
        },
        values: [],
        correction: 2,
      ),
      "gaussian": Foo(
        calc: (double n, double N, double a) {
          return pow(e, -0.5 * pow((n - (N - 1) / 2) / (a * (N - 1) / 2), 2))
              .toDouble();
        },
        values: [],
        correction: 5 / 3,
      ),
      "bartlettHanning": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.62 - 0.48 * ((n / (N - 1)) - 0.5).abs() - 0.38 * cos(z);
        },
        values: [],
        correction: 2,
      ),
      "blackman": Foo(
        calc: (double n, double N, double a) {
          var a0 = (1 - a) / 2;
          var a1 = 0.5;
          var a2 = a / 2;
          var z = (PI2 * n) / (N - 1);
          return a0 - a1 * cos(z) + a2 * cos(2 * z);
        },
        values: [],
        correction: 4 / 3,
      ),
      "blackmanHarris": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.35875 -
              0.48829 * cos(z) +
              0.14128 * cos(2 * z) -
              0.01168 * cos(3 * z);
        },
        values: [],
        correction: 1.5594508635,
      ),
      "nuttall3": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.375 - 0.5 * cos(z) + 0.125 * cos(2 * z);
        },
        values: [],
        correction: 1.56,
      ),
      "nuttall3a": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.40897 - 0.5 * cos(z) + 0.09103 * cos(2 * z);
        },
        values: [],
        correction: 1.692,
      ),
      "nuttall3b": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.4243801 - 0.4973406 * cos(z) + 0.078793 * cos(2 * z);
        },
        values: [],
        correction: 1.7372527,
      ),
      "nuttall4": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.3125 -
              0.46875 * cos(z) +
              0.1875 * cos(2 * z) -
              0.03125 * cos(3 * z);
        },
        values: [],
        correction: 1.454543,
      ),
      "nuttall4a": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.338946 -
              0.481973 * cos(z) +
              0.161054 * cos(2 * z) -
              0.018027 * cos(3 * z);
        },
        values: [],
        correction: 1.512732763,
      ),
      "nuttall4b": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.355768 -
              0.481973 * cos(z) +
              0.144232 * cos(2 * z) -
              0.012604 * cos(3 * z);
        },
        values: [],
        correction: 1.55223262,
      ),
      "nuttall4c": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.3635819 -
              0.4891775 * cos(z) +
              0.1365995 * cos(2 * z) -
              0.0106411 * cos(3 * z);
        },
        values: [],
        correction: 1.57129067,
      ),
      // fast decaying flat top
      "sft3f": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.26526 - 0.5 * cos(z) + 0.23474 * cos(2 * z);
        },
        values: [],
        correction: 1.3610238,
      ),
      "sft4f": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.21706 -
              0.42103 * cos(z) +
              0.28294 * cos(2 * z) -
              0.07897 * cos(3 * z);
        },
        values: [],
        correction: 1.2773573,
      ),
      "sft5f": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.1881 -
              0.36923 * cos(z) +
              0.28702 * cos(2 * z) -
              0.13077 * cos(3 * z) +
              0.02488 * cos(4 * z);
        },
        values: [],
        correction: 1.23167769,
      ),
      // minimum sidelobe flat top
      "sft3m": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.28235 - 0.52105 * cos(z) + 0.19659 * cos(2 * z);
        },
        values: [],
        correction: 1.39343451,
      ),
      "sft4m": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.241906 -
              0.460841 * cos(z) +
              0.2552381 * cos(2 * z) -
              0.041872 * cos(3 * z);
        },
        values: [],
        correction: 1.3190596,
      ),
      "sft5m": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.209671 -
              0.407331 * cos(z) +
              0.281225 * cos(2 * z) -
              0.092669 * cos(3 * z) +
              0.0091036 * cos(4 * z);
        },
        values: [],
        correction: 1.26529456464,
      ),
      "nift": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return 0.2810639 - 0.5208972 * cos(z) + 0.1980399 * cos(2 * z);
        },
        values: [],
        correction: 1.39094182,
      ),
      "hpft": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.912510941 * cos(z) +
                  1.079173272 * cos(2 * z) -
                  0.1832630879 * cos(3 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
      "srft": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.93 * cos(z) +
                  1.29 * cos(2 * z) -
                  0.388 * cos(3 * z) +
                  0.028 * cos(4 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
      "hft70": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.90796 * cos(z) +
                  1.07349 * cos(2 * z) -
                  0.18199 * cos(3 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
      "hft95": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.9383379 * cos(z) +
                  1.3045202 * cos(2 * z) -
                  0.402827 * cos(3 * z) +
                  0.0350665 * cos(4 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
      "hft90d": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.942604 * cos(z) +
                  1.340318 * cos(2 * z) -
                  0.440811 * cos(3 * z) +
                  0.043097 * cos(4 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
      "hft116d": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.9575375 * cos(z) +
                  1.4780705 * cos(2 * z) -
                  0.6367431 * cos(3 * z) +
                  0.1228389 * cos(4 * z) -
                  0.0066288 * cos(5 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
      "hft144d": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.96760033 * cos(z) +
                  1.57983607 * cos(2 * z) -
                  0.81123644 * cos(3 * z) +
                  0.22583558 * cos(4 * z) -
                  0.02773848 * cos(5 * z) +
                  0.0009036 * cos(6 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
      "hft196d": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.97441842 * cos(z) +
                  1.65409888 * cos(2 * z) -
                  0.95788186 * cos(3 * z) +
                  0.3367342 * cos(4 * z) -
                  0.06364621 * cos(5 * z) +
                  0.00521942 * cos(6 * z) -
                  0.00010599 * cos(7 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
      "hft223d": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.98298997309 * cos(z) +
                  1.75556083063 * cos(2 * z) -
                  1.19037717712 * cos(3 * z) +
                  0.56155440797 * cos(4 * z) -
                  0.17296769663 * cos(5 * z) +
                  0.03233247087 * cos(6 * z) -
                  0.00324954578 * cos(7 * z) +
                  0.00013801040 * cos(8 * z) -
                  0.00000132725 * cos(9 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
      "hft248d": Foo(
        calc: (double n, double N, double a) {
          var z = (PI2 * n) / (N - 1);
          return (1.0 -
                  1.985844164102 * cos(z) +
                  1.791176438506 * cos(2 * z) -
                  1.282075284005 * cos(3 * z) +
                  0.667777530266 * cos(4 * z) -
                  0.240160796576 * cos(5 * z) +
                  0.056656381764 * cos(6 * z) -
                  0.008134974479 * cos(7 * z) +
                  0.00062454465 * cos(8 * z) -
                  0.000019808998 * cos(9 * z) +
                  0.000000132974 * cos(10 * z)) /
              N;
        },
        values: [],
        correction: 1,
      ),
    };
  }

  // TODO: In JS this method returns either a number array, or a number
  // We're falling back to dynamic here in order to do this...
  dynamic windowFunctions(WinFunction params) {
    var foo = this.windowCalculation[params.name]!;
    if (foo.values.length != params.N) {
      if (params.n == 0) {
        foo.values.length = 0;
      }

      // TODO: I put this here to hopefully grow the array, but it might as well shrink it...
      foo.values.insert(
          params.n,
          foo.correction *
              foo.calc(params.n.toDouble(), params.N.toDouble(), params.a));
      return foo.values[params.n];
    }
    return foo.values;
  }

  ComplexVector forward(dynamic b, String window) {
    int n;
    int k2;
    int h;
    int d;
    double c;
    double s;
    int ik;
    double dx;
    double dy;
    n = this.fft.buffer.length;
    var winFunction = WinFunction(
      name: window,
      N: n,
      a: 0.5,
      n: 0,
    );
    var w = this.windowFunctions(winFunction);
    if (w is double) {
      for (var i = 0; i < n; ++i) {
        winFunction.n = i;
        this.fft.buffer[i] = b[i] * this.windowFunctions(winFunction);
      }
    } else {
      for (var i = 0; i < n; ++i) {
        this.fft.buffer[i] = b[i] * w[i];
      }
    }
    for (var i = n - 1; i >= 0; i--) {
      this.fft.re[i] = this.fft.buffer[this.fft.twiddle[i]];
      this.fft.im[i] = 0.0;
    }
    for (var k = 1; k < n; k = k2) {
      h = 0;
      k2 = k + k;
      d = n ~/ k2; // TODO: Was / instead of ~/
      for (var j = 0; j < k; j++) {
        c = this.fft.cosTable[h];
        s = this.fft.sinTable[h];
        for (var i = j; i < n; i += k2) {
          ik = i + k;
          dx = s * this.fft.im[ik] + c * this.fft.re[ik];
          dy = c * this.fft.im[ik] - s * this.fft.re[ik];
          this.fft.re[ik] = this.fft.re[i] - dx;
          this.fft.re[i] += dx;
          this.fft.im[ik] = this.fft.im[i] - dy;
          this.fft.im[i] += dy;
        }
        h += d;
      }
    }
    return ComplexVector(
      this.fft.re,
      this.fft.im,
    );
  }

  List<double> inverse(ComplexVector v) {
    var i;
    var j;
    var n;
    var k;
    var k2;
    int h;
    int d;
    var c;
    var s;
    var ik;
    var dx;
    var dy;
    n = v.re.length;
    for (i = n - 1; i >= 0; i--) {
      j = this.fft.twiddle[i];
      this.fft.reI[i] = v.re[j];
      this.fft.imI[i] = -v.im[j];
    }
    for (k = 1; k < n; k = k2) {
      h = 0;
      k2 = k + k;
      d = n ~/ k2;
      for (j = 0; j < k; j++) {
        c = this.fft.cosTable[h];
        s = this.fft.sinTable[h];
        for (i = j; i < n; i += k2) {
          ik = i + k;
          dx = s * this.fft.imI[ik] + c * this.fft.reI[ik];
          dy = c * this.fft.imI[ik] - s * this.fft.reI[ik];
          this.fft.reI[ik] = this.fft.reI[i] - dx;
          this.fft.reI[i] += dx;
          this.fft.imI[ik] = this.fft.imI[i] - dy;
          this.fft.imI[i] += dy;
        }
        h += d;
      }
    }
    for (i = n - 1; i >= 0; i--) {
      this.fft.buffer[i] = this.fft.reI[i] / n;
    }
    return this.fft.buffer;
  }

  magnitude(dynamic params) {
    var ret = [];
    for (var cnt = 0; cnt < params.re.length; cnt++) {
      ret.add(sqrt(
          params.re[cnt] * params.re[cnt] + params.im[cnt] * params.im[cnt]));
    }
    return ret;
  }

  magToDb(dynamic b) {
    var ret = [];
    for (var cnt = 0; cnt < b.length; cnt++) {
      ret.add(20 * log(b[cnt]) * log10e);
    }
    return ret;
  }

  phase(dynamic params) {
    var ret = [];
    for (var cnt = 0; cnt < params.re.length; cnt++) {
      ret.add(atan2(params.im[cnt], params.re[cnt]));
    }
    return ret;
  }

  windows() {
    var winFuncs = [];
    this.windowCalculation.forEach((key, value) {
      winFuncs.add(key);
    });
    return winFuncs;
  }
}

class FftData {
  int length;
  List<double> buffer;
  List<double> re;
  List<double> im;
  List<double> reI;
  List<double> imI;
  List<int> twiddle;
  List<double> sinTable;
  List<double> cosTable;
  FftData({
    required this.length,
    required this.buffer,
    required this.re,
    required this.im,
    required this.reI,
    required this.imI,
    required this.twiddle,
    required this.sinTable,
    required this.cosTable,
  });
}

class WinFunction {
  String name;
  int N;
  double a;
  int n;
  WinFunction(
      {required this.name, required this.N, required this.a, required this.n});
}
